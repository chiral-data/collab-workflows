#!/usr/bin/env python3
"""Node-level tests for workflow-029 (de novo binder design).

Covers the semantics the pipeline depends on and cannot see for itself: which
residues node 01 declares contiguous, how it resolves hotspots, and what node 02
reads back out of the target it is handed. Every expected value here was measured
against the bundled 2B5I target, so a failure means a node changed, not that the
expectation was wrong.

Node 01's tests need BioPython, which is not on a stock python3. Run them inside
the node's own image, with this workflow directory mounted:

    python3 test_workflow.py

Node 02's reader tests are stdlib-only. The ContigMap check that guards node 02's
index arithmetic needs the RFdiffusion image and is a separate script,
02_rfdiffusion/check_contig_semantics.py, not part of this run.
"""

import contextlib
import importlib.util
import io
import os
import shutil
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = os.path.join(HERE, 'input_files', 'target.pdb')
NODE_01 = os.path.join(HERE, '01_validate_target', 'script.py')
NODE_02 = os.path.join(HERE, '02_rfdiffusion', 'run.py')

OUTPUT = './outputs/validated_target.pdb'


def load_module(name, path, tolerate_exit=False):
    """Import a node script by path.

    tolerate_exit is for node 02, which has no __main__ guard: its top-level code
    runs on import and stops at the GPU check. Its functions are defined above
    that check, so they are bound by the time it fails.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except (SystemExit, OSError):
        if not tolerate_exit:
            raise
    return module


node_01 = load_module('workflow_029_node_01', NODE_01)


def ca_pdb(numbers, chain='A'):
    """A PDB body of CA records only — enough for node 02's readers, which index
    exactly those."""
    lines = []
    for serial, number in enumerate(numbers, start=1):
        lines.append(f'ATOM  {serial:>5}  CA  ALA {chain}{number:>4}    '
                     f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C')
    return '\n'.join(lines) + '\n'


class NodeTestCase(unittest.TestCase):
    """A scratch working directory shaped like a node's, since the scripts read
    ./inputs and write ./outputs relative to wherever they are run."""

    def setUp(self):
        self.origin = os.getcwd()
        self.workdir = tempfile.mkdtemp(prefix='workflow-029-')
        os.makedirs(os.path.join(self.workdir, 'inputs'))
        os.makedirs(os.path.join(self.workdir, 'outputs'))
        os.chdir(self.workdir)

    def tearDown(self):
        os.chdir(self.origin)
        shutil.rmtree(self.workdir, ignore_errors=True)

    def write_input(self, name, text):
        path = os.path.join('inputs', name)
        with open(path, 'w') as f:
            f.write(text)
        return path


class Node01TestCase(NodeTestCase):
    def setUp(self):
        super().setUp()
        shutil.copy(TARGET, os.path.join('inputs', 'target.pdb'))

    def shift_target_numbering(self, offset):
        """The bundled target's chain A alone, renumbered by offset. A negative
        shift is the one case a contig cannot express, so it is what the escape
        hatch is tested against."""
        with open(TARGET) as source, open(os.path.join('inputs', 'target.pdb'), 'w') as f:
            for line in source:
                if line.startswith(('ATOM', 'HETATM')) and line[21] == 'A':
                    number = int(line[22:26]) + offset
                    f.write(f'{line[:22]}{number:>4}{line[26:]}')
            f.write('END\n')

    def run_node_01(self, expect_exit=False, **params):
        """Run node 01 against ./inputs/target.pdb; returns its stdout."""
        settings = {
            'PARAM_RECEPTOR_CHAIN': 'A',
            'PARAM_HOTSPOT_RESIDUES': '',
            'PARAM_STRIP_HETEROATOMS': 'true',
            'PARAM_RENUMBER_RESIDUES': 'false',
        }
        settings.update({f'PARAM_{key.upper()}': value for key, value in params.items()})
        restore = {key: os.environ.get(key) for key in settings}
        os.environ.update(settings)
        stdout = io.StringIO()
        code = None
        try:
            with contextlib.redirect_stdout(stdout):
                node_01.main()
        except SystemExit as exit_signal:
            code = exit_signal.code
        finally:
            for key, value in restore.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
        if code is not None and not expect_exit:
            self.fail(f'node 01 exited {code}:\n{stdout.getvalue()}')
        if code is None and expect_exit:
            self.fail(f'node 01 was expected to exit non-zero:\n{stdout.getvalue()}')
        self.exit_code = code
        return stdout.getvalue()

    def remarks(self, keyword):
        """The REMARK 999 records of the written target carrying this keyword."""
        prefix = f'REMARK 999 {keyword}'
        with open(OUTPUT) as f:
            return [line[len(prefix):].strip() for line in f if line.startswith(prefix)]

    def segments(self):
        recorded = self.remarks('SEGMENTS')
        self.assertEqual(len(recorded), 1, 'exactly one SEGMENTS record is written')
        return recorded[0]

    def breaks(self):
        """Each BREAK as (previous, next, missing, distance), the fields the
        contig's correctness turns on; the kind is asserted where it is the point."""
        found = []
        for record in self.remarks('BREAK'):
            previous, following, missing, distance, _ = record.split()
            found.append((int(previous), int(following), int(missing), distance))
        return found

    def kinds(self):
        return [record.split()[4] for record in self.remarks('BREAK')]


class TestNode01Segments(Node01TestCase):
    """The four chains of the bundled target, as measured during the review that
    established these numbers. Node 01 relabels every chain to A on output, so a
    segment letter is always A regardless of the input chain."""

    def test_chain_a(self):
        output = self.run_node_01(receptor_chain='A')
        self.assertEqual(self.segments(), 'A6-75/A80-98/A103-133')
        self.assertEqual(self.breaks(), [(75, 80, 4, '7.24'), (98, 103, 4, '9.28')])
        self.assertIn('Residue count: 120', output)

    def test_chain_b(self):
        output = self.run_node_01(receptor_chain='B')
        self.assertEqual(self.segments(), 'A6-24/A31-207')
        self.assertEqual(self.breaks(), [(24, 31, 6, '12.51')])
        self.assertIn('Residue count: 196', output)

    def test_chain_c(self):
        output = self.run_node_01(receptor_chain='C')
        self.assertEqual(self.segments(), 'A34-224')
        self.assertEqual(self.breaks(), [])
        self.assertIn('Residue count: 191', output)

    def test_chain_d(self):
        output = self.run_node_01(receptor_chain='D')
        self.assertEqual(self.segments(), 'A1-47/A53-64/A103-165')
        self.assertEqual(self.breaks(), [(47, 53, 5, '4.87'), (64, 103, 38, '33.87')])
        self.assertIn('Residue count: 122', output)

    def test_absent_chain_is_rejected(self):
        output = self.run_node_01(receptor_chain='Z', expect_exit=True)
        self.assertEqual(self.exit_code, 1)
        self.assertIn('Chain Z not found', output)


class TestNode01Numbering(Node01TestCase):
    def test_numbering_is_original_by_default(self):
        self.run_node_01()
        self.assertEqual(self.remarks('NUMBERING'), ['ORIGINAL'])

    def test_renumbering_collapses_the_chain_to_one_segment(self):
        self.run_node_01(renumber_residues='true')
        self.assertEqual(self.remarks('NUMBERING'), ['RENUMBERED'])
        self.assertEqual(self.segments(), 'A1-120')

    def test_renumbered_breaks_survive_as_distance_only(self):
        # The gap in numbering is gone, so nothing but the measured distance is
        # left to record that these residues are not bonded.
        self.run_node_01(renumber_residues='true')
        self.assertEqual(self.breaks(), [(70, 71, 0, '7.24'), (89, 90, 0, '9.28')])
        self.assertEqual(self.kinds(), ['DISTANCE', 'DISTANCE'])

    def test_renumbering_warns_which_breaks_it_fuses(self):
        output = self.run_node_01(renumber_residues='true')
        self.assertIn('renumbering fuses 2 chain break(s)', output)
        self.assertIn('75|80', output)
        self.assertIn('98|103', output)

    def test_non_positive_numbering_names_the_escape_hatch(self):
        self.shift_target_numbering(-10)
        output = self.run_node_01(expect_exit=True)
        self.assertEqual(self.exit_code, 1)
        self.assertIn('renumber_residues', output)

    def test_renumbering_accepts_non_positive_numbering(self):
        self.shift_target_numbering(-10)
        self.run_node_01(renumber_residues='true')
        self.assertEqual(self.segments(), 'A1-120')


class TestNode01Hotspots(Node01TestCase):
    HOTSPOTS = 'A38,A42,A45,A68'

    def test_hotspots_are_an_identity_mapping_by_default(self):
        self.run_node_01(hotspot_residues=self.HOTSPOTS)
        self.assertEqual(self.remarks('HOTSPOTS'), ['A38,A42,A45,A68'])

    def test_renumbering_moves_the_hotspots_with_the_residues(self):
        self.run_node_01(hotspot_residues=self.HOTSPOTS, renumber_residues='true')
        self.assertEqual(self.remarks('HOTSPOTS'), ['A33,A37,A40,A63'])

    def test_hotspot_naming_another_chain_is_dropped(self):
        # Only the receptor chain is kept, so a hotspot on any other chain is a
        # request this node cannot honour.
        output = self.run_node_01(receptor_chain='B', hotspot_residues=self.HOTSPOTS)
        self.assertIn('names chain A, but the receptor chain is B', output)
        self.assertIn('None of the requested hotspot residues could be mapped', output)
        self.assertEqual(self.remarks('HOTSPOTS'), [''])


class TestNode02Readers(NodeTestCase):
    """Node 02's readers, which are stdlib-only and need no image."""

    BINDER_LENGTH = 80
    RECORDED = 'A6-75/A80-98/A103-133'

    def setUp(self):
        super().setUp()
        self.node_02 = load_module('workflow_029_node_02', NODE_02, tolerate_exit=True)

    def contig(self, segments):
        # The string node 02 builds; the assertion of record, since this is what
        # RFdiffusion is actually handed.
        return f'[{segments}/0 {self.BINDER_LENGTH}-{self.BINDER_LENGTH}]'

    def test_recorded_segments_are_read_back(self):
        path = self.write_input('validated_target.pdb',
                                f'REMARK 999 NUMBERING ORIGINAL\n'
                                f'REMARK 999 SEGMENTS {self.RECORDED}\n'
                                + ca_pdb([6, 7, 80, 103]))
        segments = self.node_02.read_receptor_segments(path)
        self.assertEqual(segments, self.RECORDED)
        self.assertEqual(self.contig(segments), '[A6-75/A80-98/A103-133/0 80-80]')

    def test_absent_segment_record_reads_as_none(self):
        path = self.write_input('validated_target.pdb', ca_pdb([1, 2, 3]))
        self.assertIsNone(self.node_02.read_receptor_segments(path))

    def test_segments_are_derived_when_the_record_is_missing(self):
        # A renumbered target with its remarks stripped: contiguous from 1, so the
        # derived contig says one segment and the fused breaks are invisible here.
        path = self.write_input('validated_target.pdb', ca_pdb(range(1, 121)))
        segments = self.node_02.derive_receptor_segments(path)
        self.assertEqual(segments, 'A1-120')
        self.assertEqual(self.contig(segments), '[A1-120/0 80-80]')

    def test_derivation_splits_on_every_gap(self):
        path = self.write_input('validated_target.pdb', ca_pdb([6, 7, 8, 80, 81, 103]))
        self.assertEqual(self.node_02.derive_receptor_segments(path),
                         'A6-8/A80-81/A103-103')

    def test_hotspot_mapping_is_read_back(self):
        path = self.write_input('validated_target.pdb',
                                'REMARK 999 HOTSPOT A38 A33\n'
                                'REMARK 999 HOTSPOTS A33,A37,A40,A63\n'
                                + ca_pdb([1, 2]))
        self.assertEqual(self.node_02.read_hotspot_mapping(path), 'A33,A37,A40,A63')

    def test_absent_hotspot_mapping_reads_as_none(self):
        path = self.write_input('validated_target.pdb', ca_pdb([1, 2]))
        self.assertIsNone(self.node_02.read_hotspot_mapping(path))


class TestNode02ChainSwap(NodeTestCase):
    def setUp(self):
        super().setUp()
        self.node_02 = load_module('workflow_029_node_02', NODE_02, tolerate_exit=True)

    def test_binder_becomes_chain_a(self):
        # RFdiffusion emits the receptor as A and the binder as B, following the
        # contig order; nodes 03-05 assume the reverse.
        path = os.path.join('outputs', 'design_0.pdb')
        with open(path, 'w') as f:
            f.write(ca_pdb([1, 2], chain='A'))
            f.write('TER    1000      ALA A   2\n')
            f.write(ca_pdb([1, 2], chain='B'))
            f.write('TER    1001      ALA B   2\n')
            f.write('END\n')

        self.node_02.swap_binder_to_chain_a(path)

        with open(path) as f:
            chains = [line[21] for line in f
                      if line.startswith(('ATOM', 'TER')) and len(line) > 21]
        self.assertEqual(chains, ['B', 'B', 'B', 'A', 'A', 'A'])

    def test_other_records_are_left_alone(self):
        path = os.path.join('outputs', 'design_1.pdb')
        with open(path, 'w') as f:
            f.write('REMARK   1 UNTOUCHED A\n')
            f.write(ca_pdb([1], chain='A'))
            f.write('END\n')

        self.node_02.swap_binder_to_chain_a(path)

        with open(path) as f:
            lines = f.readlines()
        self.assertEqual(lines[0], 'REMARK   1 UNTOUCHED A\n')
        self.assertEqual(lines[1][21], 'B')
        self.assertEqual(lines[-1], 'END\n')


if __name__ == '__main__':
    unittest.main()
