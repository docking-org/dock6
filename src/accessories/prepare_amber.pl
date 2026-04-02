#!/usr/bin/perl
#
# Devleena Shivakumar, TSRI 05/24/2006.
#    Master script to process a multiple ligand mol2 file, usually
#    a DOCK ranked ligands file, and a receptor PDB file into the
#    many pieces needed by DOCK's amber_score.
# Scott Brozell, TSRI 2007-8: updated error handling, added -i and -c.

my $ignore_amberize_errors = 0 ;
my $use_existing_ligand_charges = 0 ;
while ( "-" eq substr( $ARGV[0], 0, 1 ) ) {
    if ("c" eq substr( $ARGV[0], 1 ) ) {
        $use_existing_ligand_charges = 1 ;
        shift ;
    }
    elsif ("i" eq substr( $ARGV[0], 1 ) ) {
        $ignore_amberize_errors = 1 ;
        shift ;
    }
    else {
        print "\nWarning: Invalid option or DOCK_ranked_mol2_file begins with a dash:\n";
        print "    $ARGV[0]\n";
        print "    Terminating command line option parsing.\n\n";
        last ;
    }
}
usage() if ($#ARGV != 1);


########################################################################
### SECTION:1: generates receptor files; calls amberize_receptor    ####
########################################################################

$rec_pdb_file = $ARGV[1];
if ($rec_pdb_file =~ /(\w+)\.(\w+)/)
{
    $rec_file_prefix= $1;
}
open (REC, $rec_pdb_file)
    || die "\nError, cannot open Receptor PDB file: $rec_pdb_file\n";

system( "'DOCKHOME'/bin/amberize_receptor $rec_file_prefix 1> amberize_receptor.out 2>&1")
    == 0 || die "\nError from amberize_receptor; examine amberize_receptor.out\n";
print("Coordinate and parameter files for the Receptor $rec_file_prefix generated.\n");


#############################################################################
### SECTION:2: adds the AMBER_SCORE_ID to each MOLECULE in the input mol2 ###
#############################################################################

my($mol2_file) = "";
my($global_count) = 0;
my($local_count) = 0;

chomp ($mol2_file = $ARGV[0]);
open (MOL2, $mol2_file)
    || die "\nError, cannot open ligand mol2 file: $mol2_file\n";

# assign basename for output files based on input MOL2 filename:
if ($mol2_file =~ /(\w+)\.(\w+)/)
{
    $prefix = $1;
}
else
{
    $prefix = $mol2_file."_";
}

open (OUT, ">$prefix.amber_score.mol2")
    || die "\nError, cannot open output mol2 file: $prefix.amber_score.mol2\n";

while (<MOL2>)
{
    chomp($_);

    if ($_ =~ /^(@<TRIPOS>MOLECULE)/)
    {
        if( $local_count != 0 )
        {
            print OUT "@<TRIPOS>AMBER_SCORE_ID\n";
            print OUT "$prefix.$global_count\n\n\n";
        }

        print OUT "$_\n";

        $local_count++;
        $global_count++;
    }

    if ($_ =~ /^(@<TRIPOS>MOLECULE)/ && $local_count > 1)
    {
        $local_count = 1;
        next;
    }

    if ($_ !~ /^(@<TRIPOS>MOLECULE)/ && $_ !~ /^(#######)/ && $local_count == 1 )
    {
        print OUT "$_\n" ;
        next;
    }
}
print OUT "@<TRIPOS>AMBER_SCORE_ID\n";
print OUT "$prefix.$global_count\n\n\n";
print OUT "$_\n";

print("The AMBER score tagged mol2 file $prefix.amber_score.mol2 generated.\n");


##############################################################################
### SECTION:3: splits multi-MOLECULE input mol2 into individual mol2 files ###
##############################################################################

print("Splitting the multiple Ligand mol2 file into single mol2 files.\n");
print("The single mol2 files will have the prefix: $prefix \n");
my($mol2_file) = "";
my($global_count) = 0;
my($local_count) = 0;

$mol2_file = $ARGV[0];
open (MOL2, $mol2_file)
    || die "\nError, cannot open ligand mol2 file: $mol2_file\n";

# assign basename for output files based on input MOL2 filename:

while (<MOL2>)
{
    chomp($_);
    if ($_ =~ /^(@<TRIPOS>MOLECULE)/)
    {
        $global_count++;
        open (OUT, ">$prefix.$global_count.mol2")
            || die "\nError, cannot open output mol2 file: $prefix.$global_count.mol2\n";
        print OUT "$_\n";
        $local_count++;
    }

    if ($_ =~ /^(@<TRIPOS>MOLECULE)/ && $local_count > 1)
    {
        $local_count = 1;
        next;
    }

    if ($_ !~ /^(@<TRIPOS>MOLECULE)/ && $_ !~ /^(#######)/ && $local_count == 1 )
    {
        print OUT "$_\n" ;
        next;
    }
}
print("The number of single Ligand mol2 files generated is $global_count.\n");
close OUT;


###############################################################################
### SECTION:4: generate files for ligand and complex; call amberize scripts ###
###############################################################################

print("Generating coordinate and parameter files for Ligands and Complexes.\n");
if ( $use_existing_ligand_charges )
{
    print("Using existing charges from $mol2_file\n");
}
else
{
    print("Generating AM1-BCC charges.  This may be time consuming.\n");
}

for ($lignum = 1; $lignum <= $global_count; $lignum++)
{
    $count = 1;
    $sum1 = 0;
    open(smol2, "$prefix.$lignum.mol2")
        || die "\nError, cannot open mol2 file: $prefix.$lignum.mol2\n";
    while(<smol2>)
    {
        if ( /@<TRIPOS>MOLECULE/ )
        {
            $mol_name = <smol2>;
            $mol_name =~ s/\s+$//;
        }
        @s = split ' ', $_;
        $col{$count} = $s[8];
        $sum1 += $col{$count};
        $count++;
        $sum = round ($sum1);
    }
    print("Ligand $prefix.$lignum has total charge $sum\n");
    my $charge_method = '' ;
    if ( $use_existing_ligand_charges )
    {
        $charge_method = "-s 2" ;
    }
    else
    {
        $charge_method = "-nc $sum -c bcc" ;
    }

    if ( system( "'DOCKHOME'/bin/amberize_ligand $prefix.$lignum $charge_method 1> amberize_ligand.$prefix.$lignum.out 2>&1")
        == 0 )
    {
        print("Coordinate and parameter files for the Ligand $prefix.$lignum generated.\n");
    }
    else
    {
        print "\nError from amberize_ligand; the name of the ligand is\n";
        print "    $mol_name\n";
        print "    Examine amberize_ligand.$prefix.$lignum.out\n";
        $ignore_amberize_errors || die "\n";
    }

    if ( system( "'DOCKHOME'/bin/amberize_complex $rec_file_prefix $prefix.$lignum 1> amberize_complex.$prefix.$lignum.out  2>&1")
        == 0 )
    {
        print("Coordinate and parameter files for the Complex $rec_file_prefix.$prefix.$lignum generated.\n");
    }
    else
    {
        print "\nError from amberize_complex; the name of the ligand is\n";
        print "    $mol_name\n";
        print "    Examine amberize_complex.$prefix.$lignum.out\n";
        $ignore_amberize_errors || die "\n";
    }
}

print("$0 completed.\n") ;


sub round
{
    my($sum1) = shift;
    return int($sum1 + 0.5 * ($sum1 <=> 0));
}

sub usage
{
    my $usagedoc = <<EOF;
Usage: $0 [-c] [-i] DOCK_ranked_mol2_file Receptor_PDB_file

Description
       Create Amber input files from a mol2 file containing DOCKed ligands
       and from a pdb file containing the receptor consisting of protein or DNA.
       In particular, any nucleic acids are assumed to be DNAs.
       The Amber input files are necessary for the DOCK Amber Score.

Options
       -c
              Use charges from the DOCK_ranked_mol2_file.
              The default behavior is to generate AM1-BCC charges.
              This option must occur before DOCK_ranked_mol2_file.

       -i
              Ignore errors during ligand and complex amberization.
              Errors during receptor amberization still cause abortion.
              The default behavior is to abort on all errors.
              This option must occur before DOCK_ranked_mol2_file.

Operands
       DOCK_ranked_mol2_file
              the name of the ligand(s) only mol2 file.

       Receptor_PDB_file
              the name of the receptor only pdb file.
EOF
    print STDERR $usagedoc;
    exit 1;
}

exit 0;

