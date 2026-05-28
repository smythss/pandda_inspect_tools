# -*- coding: utf-8 -*-
# plugin for COOT to inspect and model pandda.analyse results
# Python 3 / GTK3 / Coot 0.9.x (CCP4 7.1) compatible version
#
# Based on the original by Tobias Krojer, MAX IV Laboratory
# MIT License - see inspect_pandda_analyse.py for full licence text

import os
import glob
import sys
import shutil
import csv
import logging

# ---------------------------------------------------------------------------
# GTK3 via gi.repository (replaces gtk / GTK2 from the original)
# ---------------------------------------------------------------------------
try:
    import gi
    gi.require_version('Gtk', '3.0')
    from gi.repository import Gtk
except ImportError:
    import gtk as Gtk          # GTK2 fallback (shouldn't be needed with CCP4 7.1)

# Orientation constants -- some Coot/CCP4 environments ship gi bindings where
# Gtk.Orientation exists as a type but its named members are not populated.
# Fall back to the underlying integer values (0=HORIZONTAL, 1=VERTICAL).
try:
    _HORIZ = Gtk.Orientation.HORIZONTAL
    _VERT  = Gtk.Orientation.VERTICAL
except AttributeError:
    _HORIZ, _VERT = 0, 1

import coot
import sys as _sys

# Python 2/3 compatibility helpers for csv file open
if _sys.version_info[0] >= 3:
    def _csv_open_r(path):
        return open(path, newline='')
    def _csv_open_w(path):
        return open(path, 'w', newline='')
else:
    def _csv_open_r(path):
        return open(path, 'rb')
    def _csv_open_w(path):
        return open(path, 'wb')


# ---------------------------------------------------------------------------
# Coot API compatibility wrappers
# ---------------------------------------------------------------------------
# In Coot 0.9.x (Python 3) several helpers that lived in __main__ moved into
# the coot module itself.  These thin wrappers insulate the rest of the code.

def _molecule_number_list():
    try:
        return coot.molecule_number_list()
    except AttributeError:
        import __main__
        return __main__.molecule_number_list()


def _set_map_displayed(imol, show):
    """Show (1) or hide (0) a map molecule."""
    try:
        coot.set_map_displayed(imol, show)
    except AttributeError:
        import __main__
        __main__.toggle_display_map(imol, show)


def _move_molecule_here(imol):
    """Move a molecule to the current screen centre / pointer."""
    try:
        coot.move_molecule_to_screen_centre_py(imol)
    except AttributeError:
        try:
            import __main__
            __main__.move_molecule_here(imol)
        except AttributeError:
            pass


def _merge_molecules(imol_list, target_imol):
    try:
        coot.merge_molecules(imol_list, target_imol)     # Coot 0.9.x
    except AttributeError:
        coot.merge_molecules_py(imol_list, target_imol)  # Coot 0.8.x


# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

def init_logger(logfile):
    logger = logging.getLogger('pandda_inspect')
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter(
        '%(asctime)s | %(levelname)s - INSPECT | %(message)s',
        '%m-%d-%Y %H:%M:%S',
    )
    if not logger.handlers:
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(fmt)
        logger.addHandler(sh)
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


# ---------------------------------------------------------------------------
# Main inspect GUI class
# ---------------------------------------------------------------------------

class inspect_gui(object):

    def __init__(self):
        self.logger = init_logger('inspect.log')
        self.logger.info('starting new session of pandda event map inspection')

        self.index = -1
        self.Todo = []
        self.cb_list = []
        self.mol_dict = {'protein': None, 'emap': None, 'ligand': None}

        self.panddaDir = None
        self.eventCSV  = None
        self.reset_params()
        self.merged = False

        self.ligand_confidence_button_labels = [
            [0, 'unassigned'],
            [1, 'no ligand bound'],
            [2, 'unknown ligand'],
            [3, 'ambiguous density'],
            [4, 'event map only'],
            [5, '2fofc map'],
        ]

        self.selection_criteria = [
            'show all events',
            'show all events - sort by cluster size',
            'show all events - sort alphabetically',
            'show not viewed events',
            'show unassigned',
            'show no ligands bound',
            'show unknown ligands',
            'show low confidence ligands',
            'show high confidence ligands',
        ]
        self.selected_selection_criterion = None

    # ------------------------------------------------------------------
    # GUI construction
    # ------------------------------------------------------------------

    def startGUI(self):
        self.window = Gtk.Window()
        self.window.connect("delete-event", lambda w, e: Gtk.main_quit())
        self.window.set_border_width(10)
        self.window.set_default_size(400, 680)
        self.window.set_title("pandda inspect")

        self.vbox = Gtk.Box(orientation=_VERT, spacing=4)

        # ---- PanDDA folder ----
        frame = Gtk.Frame(label='PanDDA folder')
        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        btn = Gtk.Button(label="Select pandda directory")
        btn.connect("clicked", self.select_pandda_folder)
        hbox.pack_start(btn, True, True, 0)
        frame.add(hbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Event selection ----
        frame = Gtk.Frame(label='Event selection')
        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        self.select_events_combobox = Gtk.ComboBoxText()
        for c in self.selection_criteria:
            self.select_events_combobox.append_text(c)
        hbox.pack_start(self.select_events_combobox, True, True, 0)
        go_btn = Gtk.Button(label="Go")
        go_btn.connect("clicked", self.select_events)
        hbox.pack_start(go_btn, False, False, 0)
        frame.add(hbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Info grid (replaces gtk.Table) ----
        outer_frame = Gtk.Frame()
        grid = Gtk.Grid()
        grid.set_row_homogeneous(False)
        grid.set_column_homogeneous(True)
        grid.set_row_spacing(2)
        grid.set_column_spacing(2)

        self.info_labels = {}
        for row, name in enumerate(
                ['Crystal', 'Resolution', 'Rwork', 'Rfree', 'Event', 'Site', 'BDC']):
            lf = Gtk.Frame()
            lf.add(Gtk.Label(label=name))
            grid.attach(lf, 0, row, 1, 1)
            val = Gtk.Label(label='')
            vf = Gtk.Frame()
            vf.add(val)
            grid.attach(vf, 1, row, 1, 1)
            self.info_labels[name] = val

        # Convenience aliases matching old attribute names
        self.xtal_label       = self.info_labels['Crystal']
        self.resolution_label = self.info_labels['Resolution']
        self.r_work_label     = self.info_labels['Rwork']
        self.r_free_label     = self.info_labels['Rfree']
        self.event_label      = self.info_labels['Event']
        self.site_label       = self.info_labels['Site']
        self.bdc_label        = self.info_labels['BDC']

        outer_frame.add(grid)
        self.vbox.pack_start(outer_frame, False, False, 0)

        # ---- Navigator ----
        frame = Gtk.Frame(label='Navigator')
        nav_vbox = Gtk.Box(orientation=_VERT, spacing=4)

        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        for label, cb in [("<<< Event", self.previous_event),
                          ("Event >>>", self.next_event)]:
            b = Gtk.Button(label=label)
            b.connect("clicked", cb)
            hbox.pack_start(b, True, True, 0)
        nav_vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        for label, cb in [("<<< Site", self.previous_site),
                          ("Site >>>", self.next_site)]:
            b = Gtk.Button(label=label)
            b.connect("clicked", cb)
            hbox.pack_start(b, True, True, 0)
        nav_vbox.pack_start(hbox, False, False, 0)

        self.cb = Gtk.ComboBoxText()
        self.cb.connect("changed", self.select_crystal)
        nav_vbox.pack_start(self.cb, False, False, 0)

        self.crystal_progressbar = Gtk.ProgressBar()
        nav_vbox.pack_start(self.crystal_progressbar, False, False, 0)

        frame.add(nav_vbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Toggle maps ----
        frame = Gtk.Frame(label='Toggle Maps')
        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        for label, cb in [("event map",    self.toggle_emap),
                          ("Z-map",        self.toggle_zmap),
                          ("(2)fofc maps", self.toggle_x_ray_maps)]:
            b = Gtk.Button(label=label)
            b.connect("clicked", cb)
            hbox.pack_start(b, True, True, 0)
        self.toggle_average_map_button = Gtk.Button(label="average map")
        self.toggle_average_map_button.connect("clicked", self.toggle_average_map)
        hbox.pack_start(self.toggle_average_map_button, True, True, 0)
        frame.add(hbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Ligand modelling ----
        frame = Gtk.Frame(label='Ligand Modeling')
        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        for label, cb in [("Place Ligand here",  self.place_ligand_here),
                          ("Merge Ligand",        self.merge_ligand_into_protein),
                          ("Revert to unfitted",  self.reset_to_unfitted)]:
            b = Gtk.Button(label=label)
            b.connect("clicked", cb)
            hbox.pack_start(b, True, True, 0)
        frame.add(hbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Annotation radio buttons ----
        frame = Gtk.Frame(label='Annotation')
        ann_vbox = Gtk.Box(orientation=_VERT, spacing=2)
        self.ligand_confidence_button_list = []
        first_radio = None
        for item in self.ligand_confidence_button_labels:
            if first_radio is None:
                btn = Gtk.RadioButton(label=item[1])
                first_radio = btn
            else:
                btn = Gtk.RadioButton.new_with_label_from_widget(
                    first_radio, item[1])
            btn.connect("toggled", self.set_ligand_confidence, item[1])
            self.ligand_confidence_button_list.append(btn)
            ann_vbox.pack_start(btn, False, False, 0)
            btn.show()
        frame.add(ann_vbox)
        self.vbox.pack_start(frame, False, False, 0)

        # ---- Save ----
        frame = Gtk.Frame(label='Save')
        hbox = Gtk.Box(orientation=_HORIZ, spacing=4)
        self.save_next_button = Gtk.Button(label="Save Model")
        self.save_next_button.connect("clicked", self.save_next)
        hbox.pack_start(self.save_next_button, True, True, 0)
        frame.add(hbox)
        self.vbox.pack_start(frame, False, False, 0)

        self.window.add(self.vbox)
        self.window.show_all()

    # ------------------------------------------------------------------
    # CSV persistence
    # ------------------------------------------------------------------

    def save_pandda_inspect_events_csv_file(self):
        path = os.path.join(self.analysis_folder, 'pandda_inspect_events.csv')
        self.logger.info('updating {0!s}'.format(path))
        with _csv_open_w(path) as f:
            csv.writer(f).writerows(self.elist)

    # ------------------------------------------------------------------
    # Annotation callbacks
    # ------------------------------------------------------------------

    def set_ligand_confidence(self, widget, data=None):
        if widget.get_active():
            self.elist[self.index][self.ligand_confidence_index] = data
            self.save_pandda_inspect_events_csv_file()

    def save_event_as_viewed(self):
        self.elist[self.index][self.viewed_index] = 'True'
        self.save_pandda_inspect_events_csv_file()

    def set_ligand_confidence_button(self):
        found = False
        for item in self.ligand_confidence_button_labels:
            if item[1] == self.ligand_confidence:
                self.ligand_confidence_button_list[item[0]].set_active(True)
                found = True
                break
        if not found:
            self.ligand_confidence_button_list[0].set_active(True)

    # ------------------------------------------------------------------
    # Folder / CSV initialisation
    # ------------------------------------------------------------------

    def select_pandda_folder(self, widget):
        dlg = Gtk.FileChooserDialog(
            title="Select PanDDA directory",
            parent=None,
            action=Gtk.FileChooserAction.SELECT_FOLDER,
        )
        dlg.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        dlg.add_button("_Open",   Gtk.ResponseType.OK)

        response = dlg.run()
        if response != Gtk.ResponseType.OK:
            dlg.destroy()
            return

        self.panddaDir = dlg.get_filename()
        dlg.destroy()

        self.analysis_folder = ''
        for candidate in ('results', 'analyses', 'analysis'):
            p = os.path.join(self.panddaDir, candidate)
            if os.path.isdir(p):
                self.analysis_folder = p
                break

        self.eventCSV = os.path.realpath(
            os.path.join(self.analysis_folder, 'pandda_inspect_events.csv'))
        self.siteCSV = os.path.realpath(
            os.path.join(self.analysis_folder, 'pandda_inspect_sites.csv'))

        if not os.path.isfile(self.eventCSV):
            analyse_csv = self.eventCSV.replace(
                'pandda_inspect_events.csv', 'pandda_analyse_events.csv')
            if not os.path.isfile(analyse_csv):
                self.logger.error('cannot find {0!s}'.format(analyse_csv))
                return
            self.initialize_inspect_events_csv_file(analyse_csv)

        if not os.path.isfile(self.eventCSV):
            self.logger.error('cannot find {0!s}'.format(self.eventCSV))
            return

        if not os.path.isfile(self.siteCSV):
            analyse_csv = self.siteCSV.replace(
                'pandda_inspect_sites.csv', 'pandda_analyse_sites.csv')
            if not os.path.isfile(analyse_csv):
                self.logger.error('cannot find {0!s}'.format(analyse_csv))
                return
            self.initialize_inspect_sites_csv_file(analyse_csv)

        if not os.path.isfile(self.siteCSV):
            self.logger.error('cannot find {0!s}'.format(self.siteCSV))
            return

        self.parsepanddaDir()

    def make_secure_copy_of_original_csv(self, csv_file):
        backup = csv_file + '.original'
        if not os.path.isfile(backup):
            self.logger.info('backing up {0!s}'.format(csv_file))
            shutil.copy(csv_file, backup)

    def initialize_inspect_events_csv_file(self, analyse_csv):
        self.make_secure_copy_of_original_csv(analyse_csv)
        with _csv_open_r(analyse_csv) as f:
            rows = list(csv.reader(f))
        for i, row in enumerate(rows):
            if i == 0:
                rows[i].extend(['Interesting', 'Ligand Placed',
                                 'Ligand Confidence', 'Comment', 'Viewed'])
            else:
                rows[i].extend(['False', 'False', 'Low', 'None', 'False'])
        out = os.path.join(self.analysis_folder, 'pandda_inspect_events.csv')
        with _csv_open_w(out) as f:
            csv.writer(f).writerows(rows)

    def initialize_inspect_sites_csv_file(self, analyse_csv):
        self.make_secure_copy_of_original_csv(analyse_csv)
        with _csv_open_r(analyse_csv) as f:
            rows = list(csv.reader(f))
        for i, row in enumerate(rows):
            if i == 0:
                rows[i].extend(['Name', 'Comment'])
            else:
                rows[i].extend(['None', 'None'])
        out = os.path.join(self.analysis_folder, 'pandda_inspect_sites.csv')
        with _csv_open_w(out) as f:
            csv.writer(f).writerows(rows)

    def parsepanddaDir(self):
        self.logger.info("reading {0!s}".format(self.eventCSV))
        with _csv_open_r(self.eventCSV) as f:
            self.elist = list(csv.reader(f))

        self.logger.info("reading {0!s}".format(self.siteCSV))
        with _csv_open_r(self.siteCSV) as f:
            self.slist = list(csv.reader(f))

        for n, item in enumerate(self.elist[0]):
            if item == 'dtag':
                self.xtal_index = n
            elif item == 'Ligand Confidence':
                self.ligand_confidence_index = n
            elif item in ('event_num', 'event_idx'):
                self.event_index = n
            elif item in ('site_num', 'site_idx'):
                self.site_index = n
            elif item in ('bdc', '1-BDC'):
                self.bdc_index = n
            elif item == 'x':
                self.x_index = n
            elif item == 'y':
                self.y_index = n
            elif item == 'z':
                self.z_index = n
            elif item in ('analysed_resolution', 'high_resolution'):
                self.resolution_index = n
            elif item == 'r_work':
                self.r_work_index = n
            elif item == 'r_free':
                self.r_free_index = n
            elif item == 'Viewed':
                self.viewed_index = n
            elif item == 'cluster_size':
                self.cluster_size_index = n
            elif item == 'Ligand Placed':
                self.ligand_placed_index = n

        self.show_content_of_event_csv_file()

    def show_content_of_event_csv_file(self):
        self.logger.info("contents of {0!s}:".format(self.eventCSV))
        for n in range(1, len(self.elist)):
            x = round(float(self.elist[n][self.x_index]), 1)
            y = round(float(self.elist[n][self.y_index]), 1)
            z = round(float(self.elist[n][self.z_index]), 1)
            self.logger.info(
                ' xtal: {0!s} - event/site: {1!s}/{2!s}'
                ' - BDC: {3!s} - x,y,z: {4!s},{5!s},{6!s}'
                ' - Res: {7!s} - Rwork/Rfree: {8!s}/{9!s}'
                ' - viewed: {10!s} - confidence: {11!s}'.format(
                    self.elist[n][self.xtal_index],
                    self.elist[n][self.event_index],
                    self.elist[n][self.site_index],
                    self.elist[n][self.bdc_index],
                    x, y, z,
                    self.elist[n][self.resolution_index],
                    self.elist[n][self.r_work_index],
                    self.elist[n][self.r_free_index],
                    self.elist[n][self.viewed_index],
                    self.elist[n][self.ligand_confidence_index],
                ))
        self.init_crystal_selection_combobox()

    # ------------------------------------------------------------------
    # Crystal / event combobox
    # ------------------------------------------------------------------

    def init_crystal_selection_combobox(self):
        self.logger.info('rebuilding crystal selection combobox')
        model = self.cb.get_model()
        if model is not None:
            while len(model) > 0:
                self.cb.remove(0)
        self.cb_list = []
        for n in range(1, len(self.elist)):
            text = '{0!s} - event: {1!s} - site: {2!s}'.format(
                self.elist[n][self.xtal_index],
                self.elist[n][self.event_index],
                self.elist[n][self.site_index],
            )
            self.cb_list.append(text)
            self.cb.append_text(text)

    def update_crystal_selection_combobox(self):
        x = self.elist[self.index][self.xtal_index]
        e = self.elist[self.index][self.event_index]
        s = self.elist[self.index][self.site_index]
        text = '{0!s} - event: {1!s} - site: {2!s}'.format(x, e, s)
        for n, i in enumerate(self.cb_list):
            if i == text:
                self.cb.set_active(n)
                break

    def select_crystal(self, widget):
        tmp = str(widget.get_active_text())
        if not tmp:
            return
        self.logger.info('selected: {0!s}'.format(tmp))
        tmpx = tmp.replace(' - event: ', ' ').replace(' - site: ', ' ')
        parts = tmpx.split()
        xtal, event, site = parts[0], parts[1], parts[2]
        index_increment = 0
        for n in range(len(self.elist)):
            if (self.elist[n][self.xtal_index] == xtal and
                    self.elist[n][self.event_index] == event and
                    self.elist[n][self.site_index] == site):
                index_increment = n - self.index
                break
        self.change_event(index_increment)

    # ------------------------------------------------------------------
    # File location helpers
    # ------------------------------------------------------------------

    def get_pdb(self, missing_files):
        pdb = ''
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        modelled = os.path.join(
            ds, 'modelled_structures',
            '{0!s}-pandda-model.pdb'.format(self.xtal))
        input_pdb = os.path.join(
            ds, '{0!s}-pandda-input.pdb'.format(self.xtal))
        if os.path.isfile(modelled):
            pdb = modelled
            self.logger.info('found pdb (modelled): {0!s}'.format(pdb))
        elif os.path.isfile(input_pdb):
            pdb = input_pdb
            self.logger.info('found pdb: {0!s}'.format(pdb))
        else:
            self.logger.error('did not find pdb file')
            missing_files = True
        return pdb, missing_files

    def load_pdb(self):
        coot.set_nomenclature_errors_on_read("ignore")
        imol = coot.handle_read_draw_molecule_with_recentre(self.pdb, 0)
        self.mol_dict['protein'] = imol
        coot.set_show_symmetry_master(1)
        coot.set_show_symmetry_molecule(imol, 1)

    def get_emap(self, missing_files):
        emap = ''
        new_pandda_output = False
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        event_number = (3 - len(str(self.event))) * '0' + str(self.event)
        candidates = [
            (os.path.join(ds, '{0!s}-event_{1!s}_1-BDC_{2!s}_map.native.mtz'.format(
                self.xtal, self.event, self.bdc)), False),
            (os.path.join(ds, '{0!s}-event_{1!s}_1-BDC_{2!s}_map.native.ccp4'.format(
                self.xtal, self.event, self.bdc)), False),
            (os.path.join(ds, '{0!s}-pandda-output-event-{1!s}.mtz'.format(
                self.xtal, event_number)), True),
        ]
        for path, is_new in candidates:
            if os.path.isfile(path):
                emap = path
                new_pandda_output = is_new
                self.logger.info('found event map: {0!s}'.format(emap))
                break
        if not emap:
            self.logger.error(
                'cannot find event map for {0!s} event {1!s}'.format(
                    self.xtal, self.event))
            missing_files = True
        return emap, new_pandda_output, missing_files

    def load_emap(self):
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.emap, "FEVENT", "PHEVENT", "1", 0, 0)
        elif self.emap.endswith(".ccp4"):
            imol = coot.read_ccp4_map(self.emap, 0)
        else:
            imol = coot.make_and_draw_map(
                self.emap, "FWT", "PHWT", "1", 0, 0)
        self.mol_dict['emap'] = imol
        coot.set_colour_map_rotation_on_read_pdb(0)
        coot.set_last_map_colour(0, 0, 1)
        self.show_emap = 1
        coot.set_contour_level_in_sigma(
            self.mol_dict['emap'], 1.0 - float(self.bdc))

    def get_zmap(self, missing_files):
        zmap = ''
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        for fname in (
            '{0!s}-z_map.native.mtz'.format(self.xtal),
            '{0!s}-z_map.native.ccp4'.format(self.xtal),
            '{0!s}-pandda-output.mtz'.format(self.xtal),
        ):
            p = os.path.join(ds, fname)
            if os.path.isfile(p):
                zmap = p
                self.logger.info('found z-map: {0!s}'.format(zmap))
                break
        if not zmap:
            self.logger.error('cannot find z-map')
            missing_files = True
        return zmap, missing_files

    def load_zmap(self):
        coot.set_default_initial_contour_level_for_difference_map(3)
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.zmap, "FZVALUES", "PHZVALUES", "1", 0, 1)
            self.mol_dict['zmap'] = imol
            coot.set_map_is_difference_map(imol, True)
        elif self.zmap.endswith(".ccp4"):
            self.mol_dict['zmap'] = coot.read_ccp4_map(self.zmap, 1)
        else:
            imol = coot.auto_read_make_and_draw_maps(self.zmap)
            self.mol_dict['zmap'] = (
                imol[0] if isinstance(imol, (list, tuple)) else imol)
            coot.set_contour_level_in_sigma(self.mol_dict['zmap'], 3)
        self.show_zmap = 1

    def get_xraymap(self, missing_files):
        p = os.path.join(
            self.panddaDir, 'processed_datasets', self.xtal,
            '{0!s}-pandda-input.mtz'.format(self.xtal))
        if os.path.isfile(p):
            self.logger.info('found xray map: {0!s}'.format(p))
            return p, missing_files
        self.logger.error('did not find xray map')
        return '', True

    def load_xraymap(self):
        imol = coot.auto_read_make_and_draw_maps(self.xraymap)
        self.mol_dict['xraymap'] = imol
        coot.set_colour_map_rotation_on_read_pdb(0)
        if isinstance(imol, (list, tuple)):
            _set_map_displayed(imol[0], self.show_xraymap)
            _set_map_displayed(imol[1], self.show_xraymap)
        else:
            _set_map_displayed(imol, self.show_xraymap)

    def get_averagemap(self):
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        for fname in (
            '{0!s}-ground-state-average-map.native.mtz'.format(self.xtal),
            '{0!s}-ground-state-average-map.native.ccp4'.format(self.xtal),
            '{0!s}-pandda-output.mtz'.format(self.xtal),
        ):
            p = os.path.join(ds, fname)
            if os.path.isfile(p):
                self.logger.info('found average map: {0!s}'.format(p))
                self.toggle_average_map_button.set_sensitive(True)
                return p
        self.logger.warning('did not find average map; disabling button')
        self.toggle_average_map_button.set_sensitive(False)
        return ''

    def load_averagemap(self):
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.zmap, "FGROUND", "PHGROUND", "1", 0, 0)
            self.mol_dict['averagemap'] = imol
        elif self.averagemap.endswith(".ccp4"):
            self.mol_dict['averagemap'] = coot.read_ccp4_map(
                self.averagemap, 0)
        else:
            imol = coot.auto_read_make_and_draw_maps(self.averagemap)
            self.mol_dict['averagemap'] = (
                imol[0] if isinstance(imol, (list, tuple)) else imol)
        coot.set_colour_map_rotation_on_read_pdb(0)
        _set_map_displayed(self.mol_dict['averagemap'], self.show_averagemap)
        coot.set_last_map_colour(0, 0, 1)

    def get_ligcif(self):
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        if self.event:
            rhofit = os.path.join(
                ds, str(self.event), 'rhofit', 'best.cif')
            if os.path.isfile(rhofit):
                self.logger.info(
                    'found ligand cif (rhofit): {0!s}'.format(rhofit))
                return rhofit
        for ligcif in glob.glob(os.path.join(ds, 'ligand_files', '*.cif')):
            self.logger.info('found ligand cif: {0!s}'.format(ligcif))
            return ligcif
        self.logger.warning('no ligand cif in: {0!s}'.format(
            os.path.join(ds, 'ligand_files')))
        return ''

    def load_ligcif(self):
        if self.ligcif and os.path.isfile(self.ligcif):
            coot.read_cif_dictionary(self.ligcif)
            pdb = self.ligcif.replace('.cif', '.pdb')
            imol = coot.handle_read_draw_molecule_with_recentre(pdb, 0)
            self.mol_dict['ligand'] = imol
            coot.set_b_factor_residue_range(imol, "X", 1, 1, 20.0)
            coot.set_occupancy_residue_range(imol, "X", 1, 1, float(self.bdc))

    # ------------------------------------------------------------------
    # Event state
    # ------------------------------------------------------------------

    def reset_params(self):
        self.xtal = self.event = self.bdc = self.site = None
        self.pdb = self.emap = self.zmap = self.xraymap = self.averagemap = None
        self.new_pandda_output = False
        self.ligcif = ''
        self.x = self.y = self.z = None
        self.resolution = self.r_free = self.r_work = None
        self.ligand_confidence = None
        self.merged = False

    def update_params(self):
        missing_files = False
        self.xtal  = self.elist[self.index][self.xtal_index]
        self.event = self.elist[self.index][self.event_index]
        self.bdc   = self.elist[self.index][self.bdc_index]
        self.site  = self.elist[self.index][self.site_index]
        self.logger.info(
            'checking files for {0!s} event {1!s} site {2!s}'.format(
                self.xtal, self.event, self.site))
        self.pdb, missing_files = self.get_pdb(missing_files)
        self.emap, self.new_pandda_output, missing_files = \
            self.get_emap(missing_files)
        self.zmap, missing_files = self.get_zmap(missing_files)
        self.xraymap, missing_files = self.get_xraymap(missing_files)
        self.averagemap = self.get_averagemap()
        self.ligcif = self.get_ligcif()
        self.x = float(self.elist[self.index][self.x_index])
        self.y = float(self.elist[self.index][self.y_index])
        self.z = float(self.elist[self.index][self.z_index])
        self.logger.info('event coords: {0!s}, {1!s}, {2!s}'.format(
            self.x, self.y, self.z))
        self.resolution = self.elist[self.index][self.resolution_index]
        self.r_free     = self.elist[self.index][self.r_free_index]
        self.r_work     = self.elist[self.index][self.r_work_index]
        self.ligand_confidence = \
            self.elist[self.index][self.ligand_confidence_index]
        return missing_files

    def update_labels(self):
        self.xtal_label.set_label(str(self.xtal))
        self.resolution_label.set_label(str(self.resolution))
        self.r_free_label.set_label(str(self.r_free))
        self.r_work_label.set_label(str(self.r_work))
        self.event_label.set_label(str(self.event))
        self.site_label.set_label(str(self.site))
        self.bdc_label.set_label(str(self.bdc))

    def recentre_on_event(self):
        coot.set_rotation_centre(self.x, self.y, self.z)

    def current_sample_matches_selection_criteria(self):
        sc = self.selected_selection_criterion
        if sc is None or sc.startswith("show all events"):
            return True
        lc = self.ligand_confidence
        if sc == "show no ligands bound":
            return "no ligand bound" in lc
        if sc == "show unknown ligands":
            return "unknown ligand" in lc
        if sc == "show low confidence ligands":
            return "low confidence" in lc
        if sc == "show high confidence ligands":
            return "high confidence" in lc
        if sc == 'show not viewed events':
            return 'True' not in self.elist[self.index][self.viewed_index]
        return False

    # ------------------------------------------------------------------
    # Refresh (loads all files for the current event index)
    # ------------------------------------------------------------------

    def RefreshData(self):
        self.reset_params()

        for imol in _molecule_number_list():
            coot.close_molecule(imol)

        self.mol_dict = {
            'pdb': None, 'emap': None, 'zmap': None,
            'ligand': None, 'xraymap': None, 'averagemap': None,
        }
        self.show_emap = self.show_zmap = \
            self.show_xraymap = self.show_averagemap = 0

        if self.index < 1:
            self.index = 1
        if self.index > len(self.elist) - 1:
            self.index = len(self.elist) - 1
            self.logger.warning('reached end of events!')
            return None

        missing_files = self.update_params()

        if self.current_sample_matches_selection_criteria() and not missing_files:
            self.logger.info('loading {0!s} event {1!s}'.format(
                self.xtal, self.event))
            self.set_ligand_confidence_button()
            self.update_labels()
            self.recentre_on_event()
            self.load_ligcif()
            self.load_pdb()
            self.load_emap()
            self.load_zmap()
            self.load_xraymap()
            if self.averagemap:
                self.load_averagemap()
            self.logger.info('setting event map as RSR map')
            coot.set_imol_refinement_map(self.mol_dict['emap'])

        elif self.current_sample_matches_selection_criteria() and missing_files:
            self.logger.error('essential files missing; skipping...')
            self.change_event(1)
        else:
            self.logger.warning(
                '{0!s} event {1!s} does not match selection; skipping'.format(
                    self.xtal, self.event))
            self.change_event(1)

    # ------------------------------------------------------------------
    # Ligand modelling callbacks
    # ------------------------------------------------------------------

    def place_ligand_here(self, widget):
        self.logger.info('moving ligand to screen centre')
        _move_molecule_here(self.mol_dict['ligand'])

    def merge_ligand_into_protein(self, widget):
        self.logger.info('merging ligand into protein')
        _merge_molecules([self.mol_dict['ligand']], self.mol_dict['protein'])
        coot.close_molecule(self.mol_dict['ligand'])
        self.merged = True

    def reset_to_unfitted(self, widget):
        for imol in _molecule_number_list():
            if 'pandda-model.pdb' in coot.molecule_name(imol):
                self.pdb = os.path.join(
                    self.panddaDir, 'processed_datasets', self.xtal,
                    '{0!s}-pandda-input.pdb'.format(self.xtal))
                coot.close_molecule(imol)
                self.load_pdb()
                break

    def check_if_modelled_structures_folder_exists(self):
        p = os.path.join(self.panddaDir, 'processed_datasets',
                         self.xtal, 'modelled_structures')
        if not os.path.isdir(p):
            self.logger.info('creating {0!s}'.format(p))
            os.mkdir(p)

    def save_next(self, widget):
        self.check_if_modelled_structures_folder_exists()
        base = os.path.join(self.panddaDir, 'processed_datasets',
                            self.xtal, 'modelled_structures')
        existing = sorted(glob.glob(os.path.join(base, 'fitted-v*.pdb')))
        if existing:
            nums = [int(os.path.basename(p)[8:12]) for p in existing]
            new_n = max(nums) + 1
        else:
            new_n = 1
        new = 'fitted-v{0:04d}.pdb'.format(new_n)
        coot.write_pdb_file(self.mol_dict['protein'],
                            os.path.join(base, new))
        model_link = os.path.join(
            base, '{0!s}-pandda-model.pdb'.format(self.xtal))
        if os.path.isfile(model_link):
            os.remove(model_link)
        os.chdir(base)
        os.system('/bin/cp {0!s} {1!s}-pandda-model.pdb'.format(
            new, self.xtal))
        if self.merged:
            self.elist[self.index][self.ligand_placed_index] = 'True'
            self.save_pandda_inspect_events_csv_file()

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def select_events(self, widget):
        self.selected_selection_criterion = \
            self.select_events_combobox.get_active_text()
        self.crystal_progressbar.set_fraction(0)
        header = self.elist[0]
        data   = self.elist[1:]
        sc = self.selected_selection_criterion
        if sc == "show all events - sort by cluster size":
            data = sorted(data, key=lambda x: x[self.cluster_size_index])
        elif sc == "show all events - sort alphabetically":
            self.logger.info("sorting alphabetically")
            data = sorted(data, key=lambda x: x[self.xtal_index])
        self.elist = [header] + data
        self.init_crystal_selection_combobox()
        self.logger.info("selection: {0!s}".format(sc))
        self.index = -1

    def previous_event(self, widget):
        self.change_event(-1)

    def next_event(self, widget):
        self.save_event_as_viewed()
        for n, b in enumerate(self.ligand_confidence_button_list):
            if b.get_active():
                for c in self.ligand_confidence_button_labels:
                    if c[0] == n:
                        self.elist[self.index][
                            self.ligand_confidence_index] = c[1]
                        self.logger.info(
                            "saving confidence as '{0!s}'".format(c[1]))
                        self.save_pandda_inspect_events_csv_file()
                        break
                break
        self.change_event(1)

    def previous_site(self, widget):
        self.logger.info('moving to previous site')
        self.change_site(-1)

    def next_site(self, widget):
        self.logger.info('moving to next site')
        self.change_site(1)

    def change_site(self, n):
        current_site = int(self.site)
        new_site = current_site + n
        index_increment = 0
        for i, item in enumerate(self.elist):
            if item[self.site_index] == str(new_site):
                index_increment = i - self.index
                break
        self.change_event(index_increment)

    def change_event(self, n):
        self.index += n
        total = max(len(self.elist) - 1, 1)
        self.crystal_progressbar.set_fraction(
            min(1.0, float(self.index) / float(total)))
        self.update_crystal_selection_combobox()
        self.RefreshData()

    # ------------------------------------------------------------------
    # Map toggle callbacks
    # ------------------------------------------------------------------

    def toggle_emap(self, widget):
        if self.mol_dict['emap'] is not None:
            self.show_emap = 1 - self.show_emap
            _set_map_displayed(self.mol_dict['emap'], self.show_emap)

    def toggle_zmap(self, widget):
        if self.mol_dict['zmap'] is not None:
            self.show_zmap = 1 - self.show_zmap
            _set_map_displayed(self.mol_dict['zmap'], self.show_zmap)

    def toggle_x_ray_maps(self, widget):
        if self.mol_dict['xraymap'] is not None:
            self.show_xraymap = 1 - self.show_xraymap
            xm = self.mol_dict['xraymap']
            if isinstance(xm, (list, tuple)):
                _set_map_displayed(xm[0], self.show_xraymap)
                _set_map_displayed(xm[1], self.show_xraymap)
            else:
                _set_map_displayed(xm, self.show_xraymap)

    def toggle_average_map(self, widget):
        if self.mol_dict['averagemap'] is not None:
            self.show_averagemap = 1 - self.show_averagemap
            _set_map_displayed(
                self.mol_dict['averagemap'], self.show_averagemap)

    def CANCEL(self, widget):
        self.window.destroy()


# ---------------------------------------------------------------------------
# Entry point  (load via Coot scripting window or .coot.py)
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    inspect_gui().startGUI()
