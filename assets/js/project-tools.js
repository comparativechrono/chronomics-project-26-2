/* =====================================================================
   The programs of this project's terminal.

   The workbench is built from the Bioinformatics workbench, which has more
   programs than this project needs. Their files were left out of
   assets/vendor/ (read mapping, RNA-seq, single-cell, somatic and
   metagenomic analyses, peak calling, gene finding and HMMER), and so were
   the Networks tab and (in this project, which reads proteins rather than
   genomes) the Genome browser tab and its `igv` command. This file takes their commands out of the terminal,
   so that `help`, `which` and the completion of names list only the
   programs that are here. A command that is not here says so.

   Kept: samtools, bcftools, bgzip, tabix, minimap2, seqtk, bedtools, jq,
   kalign, muscle, FastTree, needleman_wunsch, smith_waterman, lcs, wgsim,
   gffread, lastz, plink2, gawk, grep, sed, the GNU core utilities, diff,
   cmp, python, Rscript, curl and wget, and the shell's own commands.
   ===================================================================== */
(function () {
  'use strict';
  const MG = window.MG;
  const LEFT_OUT = [
    'fastp', 'bowtie2', 'bowtie2-align-s', 'bowtie2-build', 'bowtie2-build-s', 'bwa', 'bwa-mem2', 'hisat2', 'hisat2-build',
    'hisat2-inspect', 'stringtie', 'kallisto', 'mosdepth', 'featureCounts', 'featurecounts', 'salmon', 'cutadapt', 'lofreq',
    'bam-somaticsniper', 'somaticsniper', 'prodigal', 'hmmbuild', 'hmmsearch', 'hmmscan', 'hmmpress', 'hmmfetch', 'hmmalign',
    'hmmconvert', 'hmmemit', 'hmmstat', 'phmmer', 'jackhmmer', 'nhmmer', 'nhmmscan', 'bustools', 'macs3', 'freebayes',
    'freebayes-parallel', 'kraken2', 'kraken2-build', 'kraken2-inspect', 'bracken', 'bracken-build', 'est_abundance.py',
    'generate_kmer_distribution.py', 'kmer2read_distr', 'configureStrelkaSomaticWorkflow.py', 'runWorkflow.py', 'python2',
    'network', 'igv'
  ];
  MG.projectLeftOut = LEFT_OUT.slice();
  const T = MG.shellTools || {};
  const B = MG.shellBuiltins || {};
  for (const name of LEFT_OUT) {
    delete T[name];
    delete B[name]; // (network and igv are the shell's own commands: they opened the Networks and Genome browser tabs)
  }

  /* what the terminal says about a program that is not here (shell-wb.js wrote the Bioinformatics workbench's list;
     the messages that pointed to programs left out of this project are replaced) */
  const ABSENT = MG.shellUtil && MG.shellUtil.ABSENT;
  if (ABSENT) {
    const left = new RegExp('^(' + LEFT_OUT.map((n) => n.replace(/[.+]/g, '\\$&')).join('|') + ')$');
    const drop = ['featurecounts', 'strelka', 'kraken', 'prokka', 'esl-sfetch', 'macs2', 'kb'];
    const say = [
      ['bowtie', 'is not installed in this terminal. To align sequences to a genome there is minimap2:  minimap2 -a ref.fa reads.fq > out.sam  (-x sr for short reads; man minimap2).'],
      ['fastqc', 'is not installed in this terminal; seqtk does many small jobs on FASTA and FASTQ files (seqtk comp, seqtk seq, seqtk subseq).'],
      ['Mutect2', 'is not installed in this terminal. For variants there is bcftools.'],
      ['gatk', 'is not installed in this terminal. For variants there is bcftools (view, query, filter, norm, stats, annotate, isec, consensus); for alignments samtools; for regions bedtools; for the genotypes of many people plink2.'],
      ['blastn', 'is not installed in this terminal. BLAST searches run at EBI from the Databases tab (BLAST card). Here, for pairwise alignments there are needleman_wunsch and smith_waterman; for several sequences kalign and muscle; for genomic sequences minimap2 and lastz.'],
      ['cytoscape', 'is a desktop program and is not installed in this terminal.']
    ];
    for (let i = ABSENT.length - 1; i >= 0; i--) {
      const [re] = ABSENT[i];
      if (drop.some((n) => re.test(n))) ABSENT.splice(i, 1);
      else {
        const s = say.find(([n]) => re.test(n));
        if (s) ABSENT[i] = [re, s[1]];
      }
    }
    ABSENT.unshift([left, 'is part of the Bioinformatics workbench, but not of this project’s terminal, which has only the programs the project needs. Type  help  to see them.']);
  }
})();
