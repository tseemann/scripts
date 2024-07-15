#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Time::Piece;
#use Date::Simple;

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform read_gff3);

# "OLD GENBANK"
# http://bioperl.org/pipermail/bioperl-l/2001-April/005552.html
# ---------+---------+---------+---------+---------+---------+---------+---------
# 1       10        20        30        40        50        60        70       79
# LOCUS       AB000383     5423 bp    DNA   circular  VRL       05-FEB-1999
# 
# Positions  Contents
# ---------  --------
# 01-05      LOCUS
# 06-12      spaces
# 13-21      Locus name
# 22-22      space
# 23-29      Length of sequence, right-justified
# 31-32      bp
# 34-36      Blank, ss- (single-stranded), ds- (double-stranded), or
#            ms- (mixed-stranded)
# 37-42      Blank, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA),
#            mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA
# 43-52      Blank (implies linear) or circular

# "NEW GENBANK"
# http://www.agapow.net/science/bioinformatics/genbank-files
#LOCUS       BTV8_Netherlands_2006   1981 bp    DNA     linear   VRL 15-JUL-2008
#ACCESSION   AM498054
#DEFINITION  Bluetongue virus 8 complete viral segment 4
#
#The first line of the file has very strict positional requirements:
#It must start with 'LOCUS' and then space padding up to 12 characters
#The name follows (starting at position 13) and must not contain any space,
#but can contain punctuation and other strange characters.  The positional
#constraints of the header line thus limit the length of the name.
#To describe length, 'bp' or 'aa' must appear at position 42. This is
#preceded by the atual length and space.
#Thus the length is right-aligned in its position. All other fields are left-aligned.
#The sequence type (e.g. 'DNA') appears at position 48, but can be missing (blank).
#The sequence shape (e.g. 'linear', 'circular') appears at position 56, but can be missing.
#The division code (?) appears at position 65 and is compulsory.
#The date appears at position 69.
#There is an earlier version of Genbank, with different fixed positions. Seriously. 

# MY COMMENTS:
# locusstring=1..12 id=13..X spaces=X+1..Y size=Y+1..40 space=41 bp=42
#LOCUS       AB000383     5423 bp    DNA   circular  VRL       05-FEB-1999
#LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
#LOCUS       CP002497             1110245 bp    DNA     linear   PLN 14-NOV-2011
#LOCUS       BTV8_Netherlands_2006   1981 bp    DNA     linear   VRL 15-JUL-2008

# http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
# However, the 10 characters in the locus name are no longer sufficient 


#----------------------------------------------------------------
my %tag_map = qw(
  ID	locus_tag
  Alias	old_locus_tag
  Parent	_DELETE_
);

#----------------------------------------------------------------
my %type_map = qw(
  contig	source
  assembly_component	-
  tandem_repeat	repeat_region
  inverted_repeat	repeat_region
  signal_peptide	sig_peptide
  misc_feature	misc_feature
);
while (<DATA>) {
  my @x = split m/\t/;
  next unless @x==4;
  $type_map{ $x[1] } = $x[0];
}

#----------------------------------------------------------------
my(@Options, $verbose, $unflatten);
setOptions();

my $time = localtime;

my $gff = read_gff3(@ARGV);
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'genbank');

# create a dictionary of features to the sequence the are from
# otherwise suffer horribly if #sequences is big
my %features_from;
for my $f (@{$gff->{'features'}}) {
  push @{ $features_from{ $f->seq_id } }, $f;
}

for my $seqid (sort keys %{$gff->{'sequences'}}) {
#  inform "Converting sequence: $seqid";
  my $seq = $gff->{'sequences'}{$seqid};

  # issues with overflowing the "LOCUS" name, some hacks
  my $short_seqid = $seqid;
  if (length $short_seqid > 12) {
    $short_seqid =~ s/contig/ctg/;
    $short_seqid =~ s/scaffold/scf/;
    $short_seqid =~ s/0{2,}//;
  }
  if (length $short_seqid > 12) {
    print STDERR qq/WARNING: LOCUS id "$short_seqid" still looks to long.\n/;
  }

  # construct major container, this is the meta information
  my $gb = Bio::Seq::RichSeq->new(
    -id       => $short_seqid,
    -desc     => $seqid,
    -seq      => $seq->seq,
    -division => 'HTG',  # HTG sequences (high-throughput genomic sequences)
    -molecule => 'DNA',
    -dates    => uc( $time->strftime('%d-%b-%Y') ),  # EMBL requires 20-JUL-2008 format
  );
#  $gb->keywords('Victorian Bioinformatics Consortium', 'http://vicbioinformatics.com/');

  # add a required 'source' tag for each sequence in the GFF3
  my $source = Bio::SeqFeature::Generic->new(
    -primary =>'source', 
    -start => 1, 
    -end => $seq->length,
    -frame => 0,
    -tags => { },
  );
  $gb->add_SeqFeature($source);

  # convert all the features across
  my $count=0;
#  for my $f (@{$gff->{'features'}}) {
  for my $f ( @{$features_from{$seqid}} ) {
    print STDERR "\rProcessing | $seqid | $short_seqid | ", ++$count;

    # ensure using valid ftypes
    if (my $new = $type_map{ $f->primary_tag }) {
      $f->primary_tag($new);
    }
    # fix some tags up
    # FIXME: am i being stupid here??? going over tag_map instead of going through tags?
    for my $tag (keys %tag_map) {
      if ($f->has_tag($tag)) {
        my $newtag = $tag_map{$tag};
        my @values = $f->get_tag_values($tag);
	$f->remove_tag($tag);
	$f->add_tag_value($newtag, @values) unless $newtag eq '_DELETE_';
      }
    }
#    inform $f->get_tag_values('note') if $f->has_tag('note');
  }

  # add them all to the master object
#  $gb->add_SeqFeature( @{$gff->{'features'}} );
  $gb->add_SeqFeature( @{$features_from{$seqid}} );
  
  # write it to output, then move to next sequence
  $out->write_seq($gb);
  print STDERR "\n";
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"unflatten!", VAR=>\$unflatten, DEFAULT=>0, DESC=>"Unflatten GFF features via ID tag (*NOT WORKING*)"},
  );

#  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] genome.gff3 > genome.gbk\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
# http://www.sequenceontology.org/resources/mapping/FT_SO.txt

__DATA__
-	located_sequence_feature	SO:0000110	"-" is a placeholder for no key;  should be used when the need is merely to mark  region in order to comment on it or to use it in another feature's location;	A biological feature that can be attributed to a region of biological sequence.
-10_signal	minus_10_signal	SO:0000175	Pribnow box; a conserved region about 10 bp upstream of the start point of bacterial transcription units which may be involved in  binding RNA polymerase; consensus=TAtAaT [1,2,3,4];	A conserved region about 10-bp upstream of the start point of bacterial transcription units which may be involved in binding RNA polymerase; consensus=TAtAaT.
-35_signal	minus_35_signal	SO:0000176	a conserved hexamer about 35 bp upstream of the start point of bacterial transcription units; consensus=TTGACa or TGTTGACA;	A conserved hexamer about 35-bp upstream of the start point of bacterial transcription units; consensus=TTGACa or TGTTGACA.
3'UTR	three_prime_UTR	SO:0000205	region at the 3' end of a mature transcript (following  the stop codon) that is not translated into a protein;	A region at the 3' end of a mature transcript (following the stop codon) that is not translated into a protein.
3'clip	three_prime_clip	SO:0000557	3'-most region of a precursor transcript that is clipped off during processing;	3'-most region of a precursor transcript that is clipped off during processing.
5'UTR	five_prime_UTR	SO:0000204	region at the 5' end of a mature transcript (preceding  the initiation codon) that is not translated into a  protein;	A region at the 5' end of a mature transcript (preceding the initiation codon) that is not translated into a protein.
5'clip	five_prime_clip	SO:0000555	5'-most region of a precursor transcript that is clipped off during processing;	5' most region of a precursor transcript that is clipped off during processing.
CAAT_signal	CAAT_signal	SO:0000172	CAAT box; part of a conserved sequence located about 75 bp up-stream of the start point of eukaryotic transcription units which may be involved in RNA polymerase binding; consensus=GG(C or T)CAATCT [1,2].	Part of a conserved sequence located about 75-bp upstream of the start point of eukaryotic transcription units which may be involved in RNA polymerase binding; consensus=GG(C|T)CAATCT.
CDS	CDS	SO:0000316	coding sequence; sequence of nucleotides that corresponds with the sequence of amino acids in a protein (location includes stop codon);  feature includes amino acid conceptual translation.	A contiguous sequence which begins with, and includes, a start codon and ends with, and includes, a stop codon.
C_region	undefined		constant region of immunoglobulin light and heavy  chains, and T-cell receptor alpha, beta, and gamma  chains; includes one or more exons depending on the  particular chain	
D-loop	D_loop	SO:0000297	displacement loop; a region within mitochondrial DNA in which a short stretch of RNA is paired with one strand of DNA, displacing the original partner DNA strand in this region; also used to describe the displacement of a region of one strand of duplex DNA by a single stranded invader in the reaction catalyzed by RecA protein	Displacement loop; a region within mitochondrial DNA in which a short stretch of RNA is paired with one strand of DNA, displacing the original partner DNA strand in this region; also used to describe the displacement of a region of one strand of duplex DNA by a single stranded invader in the reaction catalyzed by RecA protein.
D_segment	D_gene	SO:0000458	Diversity segment of immunoglobulin heavy chain, and  T-cell receptor beta chain;  	germline genomic DNA including D-region with 5' UTR and 3' UTR, also designated as D-segment.
GC_signal	GC_rich_region	SO:0000173	GC box; a conserved GC-rich region located upstream of the start point of eukaryotic transcription units which may occur in multiple copies or in either orientation; consensus=GGGCGG;	A conserved GC-rich region located upstream of the start point of eukaryotic transcription units which may occur in multiple copies or in either orientation; consensus=GGGCGG.
J_segment	undefined		joining segment of immunoglobulin light and heavy  chains, and T-cell receptor alpha, beta, and gamma  chains;  	
LTR	long_terminal_repeat	SO:0000286	long terminal repeat, a sequence directly repeated at both ends of a defined sequence, of the sort typically found in retroviruses;	A sequence directly repeated at both ends of a defined sequence, of the sort typically found in retroviruses.
N_region	undefined		extra nucleotides inserted between rearranged  immunoglobulin segments.	
RBS	ribosome_entry_site	SO:0000139	ribosome binding site;	Region in mRNA where ribosome assembles.
STS	STS	SO:0000331	sequence tagged site; short, single-copy DNA sequence that characterizes a mapping landmark on the genome and can be detected by PCR; a region of the genome can be mapped by determining the order of a series of STSs;	Short (typically a few hundred base pairs) DNA sequence that has a single occurrence in a genome and whose location and base sequence are known.
S_region	undefined		switch region of immunoglobulin heavy chains;   involved in the rearrangement of heavy chain DNA leading  to the expression of a different immunoglobulin class  from the same B-cell;	
TATA_signal	TATA_box	SO:0000174	TATA box; Goldberg-Hogness box; a conserved AT-rich septamer found about 25 bp before the start point of each eukaryotic RNA polymerase II  transcript  unit which may be involved in positioning the enzyme  for correct  initiation; consensus=TATA(A or T)A(A or T) [1,2];	A conserved AT-rich septamer found about 25-bp before the start point of many eukaryotic RNA polymerase II transcript units; may be involved in positioning the enzyme for correct initiation; consensus=TATA(A|T)A(A|T).
V_region	undefined		variable region of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains;  codes for the variable amino terminal portion; can be composed of V_segments, D_segments, N_regions, and J_segments;	
V_segment	undefined		variable segment of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains; codes for most of the variable region (V_region) and the last few amino acids of the leader peptide;	
attenuator	attenuator	SO:0000140	1) region of DNA at which regulation of termination of transcription occurs, which controls the expression of some bacterial operons; 2) sequence segment located between the promoter and the first structural gene that causes partial termination of transcription	A sequence segment located between the promoter and a structural gene that causes partial termination of transcription.
conflict	undefined		independent determinations of the "same" sequence differ at this site or region; Or /compare=[accession-number.sequence-version]	
enhancer	enhancer	SO:0000165	a cis-acting sequence that increases the utilization of (some)  eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter;	A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter.
exon	exon	SO:0000147	region of genome that codes for portion of spliced mRNA,  rRNA and tRNA; may contain 5'UTR, all CDSs and 3' UTR; 	A region of the genome that codes for portion of spliced messenger RNA (SO:0000234); may contain 5'-untranslated region (SO:0000204), all open reading frames (SO:0000236) and 3'-untranslated region (SO:0000205).
gap	gap	SO:0000730	gap in the sequence	A gap in the sequence of known length. THe unkown bases are filled in with N's.
gene	gene	SO:0000704	region of biological interest identified as a gene  and for which a name has been assigned;	A locatable region of genomic sequence, corresponding to a unit of inheritance, which is associated with regulatory regions, transcribed regions and/or other functional sequence regions
iDNA	iDNA	SO:0000723	intervening DNA; DNA which is eliminated through any of several kinds of recombination;	Genomic sequence removed from the genome, as a normal event, by a process of recombination. 
intron	intron	SO:0000188	a segment of DNA that is transcribed, but removed from within the transcript by splicing together the sequences (exons) on either side of it;	A segment of DNA that is transcribed, but removed from within the transcript by splicing together the sequences (exons) on either side of it.
mRNA	mRNA	SO:0000234	messenger RNA; includes 5'untranslated region (5'UTR), coding sequences (CDS, exon) and 3'untranslated region (3'UTR);	messenger RNA is the intermediate molecule between DNA and protein. It  includes UTR and coding sequences. It does not contain introns. 
mat_peptide	mature_peptide	SO:0000419	mature peptide or protein coding sequence; coding sequence for the mature or final peptide or protein product following post-translational modification; the location does not include the stop codon (unlike the corresponding CDS);	The coding sequence for the mature or final peptide or protein product following post-translational modification.
misc_RNA	transcript	SO:0000673	any transcript or RNA product that cannot be defined by other RNA keys (prim_transcript, precursor_RNA, mRNA, 5'clip, 3'clip, 5'UTR, 3'UTR, exon, CDS, sig_peptide, transit_peptide, mat_peptide, intron, polyA_site, rRNA, tRNA, scRNA, and snRNA);	An RNA synthesized on a DNA or RNA template by an RNA polymerase.
misc_binding	binding_site	SO:0000409	site in nucleic acid which covalently or non-covalently binds another moiety that cannot be described by any other binding key (primer_bind or protein_bind);	A region on the surface of a molecule that may interact with another molecule.
misc_difference	sequence_difference	SO:0000413	feature sequence is different from that presented  in the entry and cannot be described by any other  Difference key (conflict, unsure, old_sequence,  variation, or modified_base);	A region where the sequences differs from that of a specified sequence.
misc_feature	region	SO:0000001	region of biological interest which cannot be described by any other feature key; a new or rare feature;	Continous sequence.
misc_recomb	recombination_feature	SO:0000298	site of any generalized, site-specific or replicative recombination event where there is a breakage and reunion of duplex DNA that cannot be described by other recombination keys or qualifiers of source key  (/insertion_seq, /transposon, /proviral);	
misc_signal	regulatory_region	SO:0005836	any region containing a signal controlling or altering gene function or expression that cannot be described by other signal keys (promoter, CAAT_signal, TATA_signal, -35_signal, -10_signal, GC_signal, RBS, polyA_signal, enhancer, attenuator, terminator, and rep_origin).	A DNA sequence that controls the expression of a gene.
misc_structure	sequence_secondary_structure	SO:0000002	any secondary or tertiary nucleotide structure or  conformation that cannot be described by other Structure keys (stem_loop and D-loop);	A folded sequence.
modified_base	modified_base_site	SO:0000305	the indicated nucleotide is a modified nucleotide and should be substituted for by the indicated molecule (given in the mod_base qualifier value)	A modified nucleotide, i.e. a nucleotide other than A, T, C. G or (in RNA) U.
old_sequence	undefined		the presented sequence revises a previous version of the sequence at this location; Or /compare=[accession-number.sequence-version]	
operon	operon	SO:0000178	region containing polycistronic transcript containing genes that encode enzymes that are  in the same metabolic pathway and regulatory sequences 	A group of contiguous genes transcribed as a single (polycistronic) mRNA from a single regulatory region. 
oriT	origin_of_transfer	SO:0000724	origin of transfer; region of a DNA molecule where  transfer is initiated during the process of conjugation  or mobilization	A region of a DNA molecule whre transfer is initiated during the process of conjugation or mobilization.
polyA_signal	polyA_signal_sequence	SO:0000551	recognition region necessary for endonuclease cleavage of an RNA transcript that is followed by polyadenylation; consensus=AATAAA [1];	The recognition sequence necessary for endonuclease cleavage of an RNA transcript that is followed by polyadenylation; consensus=AATAAA.
polyA_site	polyA_site	SO:0000553	site on an RNA transcript to which will be added adenine residues by post-transcriptional polyadenylation;	The site on an RNA transcript to which will be added adenine residues by post-transcriptional polyadenylation.
precursor_RNA	primary_transcript	SO:0000185	any RNA species that is not yet the mature RNA product; may include 5' clipped region (5'clip), 5' untranslated region (5'UTR), coding sequences (CDS, exon), intervening sequences (intron), 3' untranslated region (3'UTR), and 3' clipped region (3'clip);	The primary (initial, unprocessed) transcript; includes five_prime_clip (SO:0000555), five_prime_untranslated_region (SO:0000204), open reading frames (SO:0000236), introns (SO:0000188) and three_prime_ untranslated_region (three_prime_UTR), and three_prime_clip (SO:0000557).
prim_transcript	primary_transcript	SO:0000185	primary (initial, unprocessed) transcript;  includes 5' clipped region (5'clip), 5' untranslated region (5'UTR), coding sequences (CDS, exon), intervening sequences (intron), 3' untranslated region (3'UTR), and 3' clipped region (3'clip);	The primary (initial, unprocessed) transcript; includes five_prime_clip (SO:0000555), five_prime_untranslated_region (SO:0000204), open reading frames (SO:0000236), introns (SO:0000188) and three_prime_ untranslated_region (three_prime_UTR), and three_prime_clip (SO:0000557).
primer_bind	primer_binding_site	SO:0005850	non-covalent primer binding site for initiation of replication, transcription, or reverse transcription; includes site(s) for synthetic e.g., PCR primer elements;	Non-covalent primer binding site for initiation of replication, transcription, or reverse transcription.  
promoter	promoter	SO:0000167	region on a DNA molecule involved in RNA polymerase binding to initiate transcription;	The region on a DNA molecule involved in RNA polymerase binding to initiate transcription.
protein_bind	protein_binding_site	SO:0000410	non-covalent protein binding site on nucleic acid;	A region of a molecule that binds to a protein.
rRNA	rRNA	SO:0000252	mature ribosomal RNA ; RNA component of the ribonucleoprotein particle (ribosome) which assembles amino acids into proteins.	RNA that comprises part of a ribosome, and that can provide both structural scaffolding and catalytic activity.
repeat_region	repeat_region	SO:0000657	region of genome containing repeating units;	A region of sequence containing one or more repeat units.
repeat_unit	repeat_unit	SO:0000726	single repeat element;	A single repeat element.
satellite	satellite_DNA	SO:0000005	many tandem repeats (identical or related) of a short basic repeating unit;  many have a base composition or other property different from the genome average  that allows them to be separated from the bulk (main band) genomic DNA;	The many tandem repeats (identical or related) of a short basic repeating unit; many have a base composition or other property different from the genome average that allows them to be separated from the bulk (main band) genomic DNA.
scRNA	scRNA	SO:0000013	small cytoplasmic RNA; any one of several small cytoplasmic RNA molecules present in the cytoplasm and (sometimes) nucleus of a eukaryote; 	Any one of several small cytoplasmic RNA moleculespresent in the cytoplasm and sometimes nucleus of a eukaryote. 
sig_peptide	signal_peptide	SO:0000418	signal peptide coding sequence; coding sequence for an N-terminal domain of a secreted protein; this domain is involved in attaching nascent polypeptide to the membrane leader sequence;	The sequence for an N-terminal domain of a secreted protein; this domain is involved in attaching nascent polypeptide to the membrane leader sequence.
snRNA	snRNA	SO:0000274	small nuclear RNA molecules involved in pre-mRNA  splicing and processing  	Small non-coding RNA in the nucleoplasm.
snoRNA	snoRNA	SO:0000275	small nucleolar RNA molecules mostly involved in  rRNA modification and processing;  	Small nucleolar RNAs (snoRNAs) are involved in the processing and modification of rRNA in the nucleolus. There are two main classes of snoRNAs: the box C/D class, and the box H/ACA class. U3 snoRNA is a member of the box C/D class. Indeed, the box C/D element is a subset of the six short sequence elements found in all U3 snoRNAs, namely boxes A, A', B, C, C', and D. The U3 snoRNA secondary structure is characterised by a small 5' domain (with boxes A and A'), and a larger 3' domain (with boxes B, C, C', and D), the two domains being linked by a single-stranded hinge. Boxes B and C form the B/C motif, which appears to be exclusive to U3 snoRNAs, and boxes C' and D form the C'/D motif. The latter is functionally similar to the C/D motifs found in other snoRNAs. The 5' domain and the hinge region act as a pre-rRNA-binding domain. The 3' domain has conserved protein-binding sites. Both the box B/C and box C'/D motifs are sufficient for nuclear retention of U3 snoRNA. The box C'/D motif is also necessary for nucleolar localization, stability and hypermethylation of U3 snoRNA. Both box B/C and C'/D motifs are involved in specific protein interactions and are necessary for the rRNA processing functions of U3 snoRNA.
source	databank_entry	SO:2000061	identifies the biological source of the specified span of the sequence; this key is mandatory; more than one source key per sequence is allowed; every entry/record will have, as a minimum, either a single source key spanning the  entire sequence or multiple source keys which together  span the entire sequence. /mol_type="genomic DNA", "genomic RNA", "mRNA",  "tRNA", "rRNA", "snoRNA", "snRNA", "scRNA", "pre-RNA", "other RNA", "other DNA",  "unassigned DNA", "unassigned RNA"	The sequence referred to by an entry in a databank such as Genbank or SwissProt.
stem_loop	stem_loop	SO:0000313	hairpin; a double-helical region formed by base-pairing between adjacent (inverted) complementary sequences in a single strand of RNA or DNA. 	A double-helical region of nucleic acid formed by base-pairing between adjacent (inverted) complementary sequences.
tRNA	tRNA	SO:0000253	mature transfer RNA, a small RNA molecule (75-85 bases long) that mediates the translation of a nucleic acid sequence into an amino acid sequence;	Transfer RNA (tRNA) molecules are approximately 80 nucleotides in length. Their secondary structure includes four short double-helical elements and three loops (D, anti-codon, and T loops). Further hydrogen bonds mediate the characteristic L-shaped molecular structure. tRNAs have two regions of fundamental functional importance: the anti-codon, which is responsible for specific mRNA codon recognition, and the 3' end, to which the tRNA's corresponding amino acid is attached (by aminoacyl-tRNA synthetases). tRNAs cope with the degeneracy of the genetic code in two manners: having more than one tRNA (with a specific anti-codon) for a particular amino acid; and 'wobble' base-pairing, i.e. permitting non-standard base-pairing at the 3rd anti-codon position.
terminator	terminator	SO:0000141	sequence of DNA located either at the end of the transcript that causes RNA polymerase to terminate  transcription;	The sequence of DNA located either at the end of the transcript that causes RNA polymerase to terminate transcription.
transit_peptide	transit_peptide	SO:0000725	transit peptide coding sequence; coding sequence for an N-terminal domain of a nuclear-encoded organellar protein; this domain is involved in post-translational import of the protein into the organelle;	The coding sequence for an N-terminal domain of a nuclear-encoded organellar protein: this domain is involved in post translational import of the protein into the organelle. 
unsure	undefined		author is unsure of exact sequence in this region;	
variation	sequence_variant	SO:0000109	a related strain contains stable mutations from the same gene (e.g., RFLPs, polymorphisms, etc.) which differ from the presented sequence at this location (and possibly others);	A region of sequence where variation has been observed.
