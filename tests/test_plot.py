
import tempfile
import unittest
import os
import hashlib

from PIL import Image

from shabam import seqplot

def checksum_file(path):
    with open(path, 'rb') as handle:
        return hashlib.sha1(handle.read()).hexdigest()

class TestPlot(unittest.TestCase):
    
    def test_seqplot(self):
        ''' test that we can make plots
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.bam')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        # check that we make a plot with the correct image. check this by
        # returning png data, then getting the sha1 digest. This will fail if
        # a single bit chenges, so perhaps this is too stringent.
        first = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=first.name)
        
        # test we don't get the empty file checksum, the file has some content
        self.assertNotEqual(checksum_file(first.name), 'da39a3ee5e6b4b0d3255bfef95601890afd80709')
        
        # now try writing a plot to PNG file.
        single = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=single.name)
        
        self.assertTrue(os.path.getsize(single.name) > 0)
        # check we have written the same data between runs
        self.assertEqual(checksum_file(first.name), checksum_file(single.name))
        
        # now check we can plot multiple sequence files in a single image
        double = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path, path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=double.name)
        self.assertTrue(os.path.getsize(double.name) > os.path.getsize(single.name))
        
        # check that the sha1 hash does not match the earlier one. this may
        # be excessive.
        self.assertNotEqual(checksum_file(single.name), checksum_file(double.name))
    
    def test_seqplot_png(self):
        ''' test that a seqplot PNG is the right shape
        '''
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.bam')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        first = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=first.name)
        
        # test basic characteristics of the written png
        im = Image.open(first.name)
        self.assertEqual(im.size, (4000, 615))
    
    def test_seqplot_by_strand(self):
        ''' test that we can plot bam, and color reads by strand
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.bam')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        # check that we can plot by strand
        data = seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            by_strand=True)
        checksum = hashlib.sha1(data).hexdigest()
        
        data = seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta)
        checksum_unstranded = hashlib.sha1(data).hexdigest()
        self.assertNotEqual(checksum, checksum_unstranded)
    
    def test_seqplot_cram(self):
        ''' test that we can plot from CRAM files
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.cram')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        single = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=single.name)
        
        # check that we have written some data to the file
        self.assertTrue(os.path.getsize(single.name) > 0)
