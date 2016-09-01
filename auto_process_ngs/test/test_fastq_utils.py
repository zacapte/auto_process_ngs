#######################################################################
# Tests for fastq_utils.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from bcftbx.mock import MockIlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from auto_process_ngs.fastq_utils import assign_barcodes_single_end
from auto_process_ngs.fastq_utils import collect_fastqs

fastq_r1 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
"""

fastq_r1_out = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TTTAC
AACTAGCTTCTCTTTTTCTT
+
31@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TCCCC
AGTCTCAGCCCTACTCCACT
+
11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:CCACC
ACGCCTGGCTAATTTTTTTT
+
AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:CAGCA
ATATACACTTCACTCTGCAT
+
1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:GATAA
AGACAGAGTCTTAATTAAAC
+
11DFFCFFDGGGB3BF313A
"""

# assign_barcodes_single_end
class TestAssignBarcodesSingleEnd(unittest.TestCase):
    """Tests for the assign_barcodes_single_end function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test file
        self.fastq_in = os.path.join(self.wd,'test.fq')
        open(self.fastq_in,'w').write(fastq_r1)
        # Output file
        self.fastq_out = os.path.join(self.wd,'out.fq')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_assign_barcodes_single_end(self):
        """assign_barcodes_single_end: extract inline barcodes from first 5 bases
        """
        nreads = assign_barcodes_single_end(self.fastq_in,
                                            self.fastq_out)
        self.assertEqual(nreads,5)
        self.assertEqual(open(self.fastq_out,'r').read(),
                         fastq_r1_out)

# collect_fastqs
class TestCollectFastqs(unittest.TestCase):
    """Tests for the collect_fastqs function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_collect_fastqs')

    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_collect_from_casava_paired_end(self):
        """collect_fastqs: test with paired-end CASAVA output
        """
        # Make an example paired-end CASAVA output dir
        mockdata = MockIlluminaData('160824_M0001_00113_ABC123','casava',
                                    unaligned_dir='casava',
                                    paired_end=True,
                                    top_dir=self.wd)
        mockdata.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=(1,2))
        mockdata.add_fastq_batch('AB','AB2','AB2_GCCAAT',lanes=(1,2))
        mockdata.add_fastq_batch('CD','CD3','CD3_ACGCCA',lanes=(3,4))
        mockdata.add_fastq_batch('CD','CD4','CD4_ACGCCA',lanes=(3,4))
        mockdata.add_undetermined(lanes=(1,2,3,4))
        mockdata.create()
        # Collect all fastqs
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in (
                'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R2_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L001_R2_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R2_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L003_R2_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L004_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L004_R2_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L003_R2_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L004_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L004_R2_001.fastq.gz',
                'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R2_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz',
                'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R2_001.fastq.gz',
                'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R2_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect all fastqs directly from 'unaligned' dir
        fastqs = collect_fastqs(mockdata.unaligned_dir)
        self.assertEqual(fastqs,expected_fastqs)
        # Collect fastqs from lane 2
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in (
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R2_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                lane=2,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect R1 fastqs only
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in (
                'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L004_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L004_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                read=1,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)

    def test_collect_from_casava_single_end(self):
        """collect_fastqs: test with single end CASAVA output
        """
        # Make an example single-end CASAVA output dir
        mockdata = MockIlluminaData('160824_M0001_00113_ABC123','casava',
                                    unaligned_dir='casava',
                                    paired_end=False,
                                    top_dir=self.wd)
        mockdata.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=(1,2))
        mockdata.add_fastq_batch('AB','AB2','AB2_GCCAAT',lanes=(1,2))
        mockdata.add_fastq_batch('CD','CD3','CD3_ACGCCA',lanes=(3,4))
        mockdata.add_fastq_batch('CD','CD4','CD4_ACGCCA',lanes=(3,4))
        mockdata.add_undetermined(lanes=(1,2,3,4))
        mockdata.create()
        # Collect all fastqs
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in (
                'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L001_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD3/CD3_ACGCCA_L004_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L003_R1_001.fastq.gz',
                'Project_CD/Sample_CD4/CD4_ACGCCA_L004_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect R1 fastqs only
        fastqs = collect_fastqs(mockdata.dirn,
                                read=1,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect R2 fastqs only
        fastqs = collect_fastqs(mockdata.dirn,
                                read=2,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,[])
        # Collect all fastqs directly from 'unaligned' dir
        fastqs = collect_fastqs(mockdata.unaligned_dir)
        self.assertEqual(fastqs,expected_fastqs)
        # Collect fastqs from lane 2
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in (
                'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                'Project_AB/Sample_AB2/AB2_GCCAAT_L002_R1_001.fastq.gz',
                'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                lane=2,
                                unaligned_dir='casava')
        self.assertEqual(fastqs,expected_fastqs)

    def test_collect_from_bcl2fastq2_paired_end(self):
        """collect_fastqs: test with paired-end bcl2fastq2 output
        """
        # Make an example paired-end bcl2fastq v2 output dir
        mockdata = MockIlluminaData('160824_M0001_00113_ABC123','bcl2fastq2',
                                    unaligned_dir='bcl2fastq2',
                                    paired_end=True,
                                    top_dir=self.wd)
        mockdata.add_fastq_batch('AB','AB1','AB1_S1',lanes=(1,2))
        mockdata.add_fastq_batch('AB','AB2','AB2_S2',lanes=(1,2))
        mockdata.add_fastq_batch('CD','CD3','CD3_S3',lanes=(3,4))
        mockdata.add_fastq_batch('CD','CD4','CD4_S4',lanes=(3,4))
        mockdata.add_undetermined(lanes=(1,2,3,4))
        mockdata.create()
        # Collect all fastqs
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                                     'AB/AB1_S1_L001_R2_001.fastq.gz',
                                     'AB/AB1_S1_L002_R1_001.fastq.gz',
                                     'AB/AB1_S1_L002_R2_001.fastq.gz',
                                     'AB/AB2_S2_L001_R1_001.fastq.gz',
                                     'AB/AB2_S2_L001_R2_001.fastq.gz',
                                     'AB/AB2_S2_L002_R1_001.fastq.gz',
                                     'AB/AB2_S2_L002_R2_001.fastq.gz',
                                     'CD/CD3_S3_L003_R1_001.fastq.gz',
                                     'CD/CD3_S3_L003_R2_001.fastq.gz',
                                     'CD/CD3_S3_L004_R1_001.fastq.gz',
                                     'CD/CD3_S3_L004_R2_001.fastq.gz',
                                     'CD/CD4_S4_L003_R1_001.fastq.gz',
                                     'CD/CD4_S4_L003_R2_001.fastq.gz',
                                     'CD/CD4_S4_L004_R1_001.fastq.gz',
                                     'CD/CD4_S4_L004_R2_001.fastq.gz',
                                     'Undetermined_S0_L001_R1_001.fastq.gz',
                                     'Undetermined_S0_L001_R2_001.fastq.gz',
                                     'Undetermined_S0_L002_R1_001.fastq.gz',
                                     'Undetermined_S0_L002_R2_001.fastq.gz',
                                     'Undetermined_S0_L003_R1_001.fastq.gz',
                                     'Undetermined_S0_L003_R2_001.fastq.gz',
                                     'Undetermined_S0_L004_R1_001.fastq.gz',
                                     'Undetermined_S0_L004_R2_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect all fastqs directly from 'unaligned' dir
        fastqs = collect_fastqs(mockdata.unaligned_dir)
        self.assertEqual(fastqs,expected_fastqs)
        # Collect fastqs from lane 2
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_L002_R1_001.fastq.gz',
                                     'AB/AB1_S1_L002_R2_001.fastq.gz',
                                     'AB/AB2_S2_L002_R1_001.fastq.gz',
                                     'AB/AB2_S2_L002_R2_001.fastq.gz',
                                     'Undetermined_S0_L002_R1_001.fastq.gz',
                                     'Undetermined_S0_L002_R2_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                lane=2,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect R1 fastqs only
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                                     'AB/AB1_S1_L002_R1_001.fastq.gz',
                                     'AB/AB2_S2_L001_R1_001.fastq.gz',
                                     'AB/AB2_S2_L002_R1_001.fastq.gz',
                                     'CD/CD3_S3_L003_R1_001.fastq.gz',
                                     'CD/CD3_S3_L004_R1_001.fastq.gz',
                                     'CD/CD4_S4_L003_R1_001.fastq.gz',
                                     'CD/CD4_S4_L004_R1_001.fastq.gz',
                                     'Undetermined_S0_L001_R1_001.fastq.gz',
                                     'Undetermined_S0_L002_R1_001.fastq.gz',
                                     'Undetermined_S0_L003_R1_001.fastq.gz',
                                     'Undetermined_S0_L004_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                read=1,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)

    def test_collect_from_bcl2fastq2_single_end(self):
        """collect_fastqs: test with single end bcl2fastq2 output
        """
        # Make an example single-end bcl2fastq v2 output dir
        mockdata = MockIlluminaData('160824_M0001_00113_ABC123','bcl2fastq2',
                                    unaligned_dir='bcl2fastq2',
                                    paired_end=False,
                                    top_dir=self.wd)
        mockdata.add_fastq_batch('AB','AB1','AB1_S1',lanes=(1,2))
        mockdata.add_fastq_batch('AB','AB2','AB2_S2',lanes=(1,2))
        mockdata.add_fastq_batch('CD','CD3','CD3_S3',lanes=(3,4))
        mockdata.add_fastq_batch('CD','CD4','CD4_S4',lanes=(3,4))
        mockdata.add_undetermined(lanes=(1,2,3,4))
        mockdata.create()
        # Collect all fastqs
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                                     'AB/AB1_S1_L002_R1_001.fastq.gz',
                                     'AB/AB2_S2_L001_R1_001.fastq.gz',
                                     'AB/AB2_S2_L002_R1_001.fastq.gz',
                                     'CD/CD3_S3_L003_R1_001.fastq.gz',
                                     'CD/CD3_S3_L004_R1_001.fastq.gz',
                                     'CD/CD4_S4_L003_R1_001.fastq.gz',
                                     'CD/CD4_S4_L004_R1_001.fastq.gz',
                                     'Undetermined_S0_L001_R1_001.fastq.gz',
                                     'Undetermined_S0_L002_R1_001.fastq.gz',
                                     'Undetermined_S0_L003_R1_001.fastq.gz',
                                     'Undetermined_S0_L004_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect all fastqs directly from 'unaligned' dir
        fastqs = collect_fastqs(mockdata.unaligned_dir)
        self.assertEqual(fastqs,expected_fastqs)
        # Collect R1 fastqs only
        fastqs = collect_fastqs(mockdata.dirn,
                                read=1,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)
        fastqs = collect_fastqs(mockdata.dirn,
                                read=2,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,[])
        # Collect fastqs from lane 2
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_L002_R1_001.fastq.gz',
                                     'AB/AB2_S2_L002_R1_001.fastq.gz',
                                     'Undetermined_S0_L002_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                lane=2,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)

    def test_collect_from_bcl2fastq2_no_lane_splitting(self):
        """collect_fastqs: test with paired-end bcl2fastq2 output (no lane splitting)
        """
        # Make example paired-end bcl2fastqv2 output dir with no lane splitting
        mockdata = MockIlluminaData('160824_M0001_00113_ABC123','bcl2fastq2',
                                    unaligned_dir='bcl2fastq2',
                                    paired_end=True,
                                    no_lane_splitting=True,
                                    top_dir=self.wd)
        mockdata.add_fastq_batch('AB','AB1','AB1_S1',lanes=(1,2))
        mockdata.add_fastq_batch('AB','AB2','AB2_S2',lanes=(1,2))
        mockdata.add_fastq_batch('CD','CD3','CD3_S3',lanes=(3,4))
        mockdata.add_fastq_batch('CD','CD4','CD4_S4',lanes=(3,4))
        mockdata.add_undetermined(lanes=(1,2,3,4))
        mockdata.create()
        # Collect all fastqs
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_R1_001.fastq.gz',
                                     'AB/AB1_S1_R2_001.fastq.gz',
                                     'AB/AB2_S2_R1_001.fastq.gz',
                                     'AB/AB2_S2_R2_001.fastq.gz',
                                     'CD/CD3_S3_R1_001.fastq.gz',
                                     'CD/CD3_S3_R2_001.fastq.gz',
                                     'CD/CD4_S4_R1_001.fastq.gz',
                                     'CD/CD4_S4_R2_001.fastq.gz',
                                     'Undetermined_S0_R1_001.fastq.gz',
                                     'Undetermined_S0_R2_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)
        # Collect all fastqs directly from 'unaligned' dir
        fastqs = collect_fastqs(mockdata.unaligned_dir)
        self.assertEqual(fastqs,expected_fastqs)
        # Collect fastqs from lane 2
        fastqs = collect_fastqs(mockdata.dirn,
                                lane=2,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,[])
        # Collect R1 fastqs only
        expected_fastqs = [os.path.join(mockdata.unaligned_dir,f)
                           for f in ('AB/AB1_S1_R1_001.fastq.gz',
                                     'AB/AB2_S2_R1_001.fastq.gz',
                                     'CD/CD3_S3_R1_001.fastq.gz',
                                     'CD/CD4_S4_R1_001.fastq.gz',
                                     'Undetermined_S0_R1_001.fastq.gz')]
        fastqs = collect_fastqs(mockdata.dirn,
                                read=1,
                                unaligned_dir='bcl2fastq2')
        self.assertEqual(fastqs,expected_fastqs)

    def test_collect_from_non_bcl2fastq_dir(self):
        """collect_fastqs: test with non-CASAVA/bcl2fastq directory
        """
        self.assertRaises(IlluminaDataError,collect_fastqs,self.wd)
