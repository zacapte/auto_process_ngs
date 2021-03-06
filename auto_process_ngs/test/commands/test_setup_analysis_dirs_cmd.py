#######################################################################
# Tests for setup_analysis_dirs_cmd.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.setup_analysis_dirs_cmd import setup_analysis_dirs

# Unit tests

class TestSetupAnalysisDirs(unittest.TestCase):
    """
    Tests for the 'setup_analysis_dirs' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSetupAnalysisDirsCommand')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0755)
            os.chmod(name,0655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_setup_analysis_dirs(self):
        """
        setup_analysis_dirs: test create new analysis dirs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_missing_metadata(self):
        """
        setup_analysis_dirs: raise exception if metadata not set
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Attempt to set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        self.assertRaises(Exception,
                          setup_analysis_dirs,ap)

    def test_setup_analysis_dirs_ignore_missing_metadata(self):
        """
        setup_analysis_dirs: test create new analysis dirs (ignore missing metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap,ignore_missing_metadata=True)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
