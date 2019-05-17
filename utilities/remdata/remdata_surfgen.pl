#!/usr/bin/perl

# Switch atoms, generating symmetric data for surfgen.

use strict;
use warnings;

my $error_flag = 0;
my $usage = "remdata_surfgen.pl [natoms] [nstates] [inputfile]\n";
my $info  = "Input file contains points to remove on each line.";

# Get command line input
$error_flag = $error_flag + &check_command_line_argument_number(3, @ARGV);
if ($error_flag != 0) {
    print "Usage: $usage$info\n";
    exit;
}
my ($natoms, $nstates, $infile) = &read_command_line_arguments(4, @ARGV);

# Check current directory for necessary data files
$error_flag = $error_flag + &check_files_exist($nstates);
if ($error_flag != 0) { exit; }

# Get number of geometries and print value.
my $ngeoms = &get_geometry_number();
print "Number of starting geometries: $ngeoms\n";

# Read input file for points to remove from data set.
$error_flag = &check_for_file($infile);
if ($error_flag != 0) { exit; }
my (@rempts) = &read_lines_of_file($infile);
if ($#rempts + 1 <= 0) { print "No points to remove!"; exit; }

# Remove data points from files.
$error_flag = &remove_datapts_energy(@rempts);
if ($error_flag != 0) { exit; }
for (my $i = 1; $i <= $nstates; $i++) {
    $error_flag = &remove_datapts_grad("cartgrd.drt1.state${i}.all", 
                                       $natoms, @rempts);
    if ($error_flag != 0) { exit; }
    for (my $j = ($i + 1); $j <= $nstates; $j++) {
        $error_flag = &remove_datapts_grad(
            "cartgrd_total.drt1.state${i}.drt1.state${j}.all",
            $natoms, @rempts);
        if ($error_flag != 0) { exit; }
    }
}


#-------------------------------------------------------------------------------
# Subroutines:
#  check_command_line_argument_number
#  check_files_exist
#  get_geometry_number
#  read_command_line_arguments
#  read_lines_of_file
#  remove_datapts_energy
#  remove_datapts_grad
#-------------------------------------------------------------------------------

# check_command_line_argument_number: check number of arguments
sub check_command_line_argument_number {
    my ($num, @args) = @_;
    if (($#args + 1) != $num) {
	print "Incorrect argument number: $#args != $num\n";
	return 1;
    } else {
	return 0;
    }
}

# check_files_exist: check the the following necessary files exist:
#  geom.all, energy.all, cartgrd.drt1.state$n.all, 
#  cartgrd_total.drt1.state$n.drt1.state$m.all
sub check_files_exist {
    my ($nst) = @_;
    if (!-e "geom.all") { print "Missing file: geom.all\n"; return 1;}
    if (!-e "energy.all") { print "Missing file: energy.all\n"; return 1;}
    for (my $i = 1; $i <= $nst; $i++) {
	if (!-e "cartgrd.drt1.state${i}.all") {
	    print "Missing file: cartgrd.drt1.state${i}.all\n"; return 1;
	}
	for (my $j = ($i + 1); $j <= $nst; $j++) {
	    if (!-e "cartgrd_total.drt1.state${i}.drt1.state${j}.all") {
		print "Missing file: cartgrd_total.drt1.state${i}.drt1.state${j}.all\n";
		return 1;
	    }
	}
    }
    return 0;
}

# check_for_file: check to see if a file exists
sub check_for_file {
    my ($name) = @_;
    if (!-e "$name") { print "Missing file: $name\n"; return 1;}
    return 0;
}

# get_geometry_number: get number of geometries
sub get_geometry_number {
    open(FILE, "<energy.all");
    my @array = grep(/\n/i, <FILE>);
    close(FILE);
    return ($#array + 1);
}

# read_command_line_arguments: read command line arguments
sub read_command_line_arguments {
    my ($num, @args) = @_;
    return @args;
}

# read_lines_of_file: read lines of a file
sub read_lines_of_file {
    my ($file) = @_;
    open(FILE, "<$file") or die "Could not open file: $file!\n";
    my @vals = grep(/\n/i, <FILE>);
    chomp(@vals);
    close(FILE);
    return @vals;
}

# remove_datapts_energy: remove data points from energy file
sub remove_datapts_energy {
    my (@pts) = @_;

    open(FILE, "<energy.all") or die "Could not open energy.all!\n";
    my @data = grep(/\n/i, <FILE>);
    close(FILE);

    system("mv energy.all energy.all.old");
    open(FILE, ">energy.all");
    my $index;
    my $flag;
    my $ptnum = 0;
    for (my $i = 0; $i <= $#data; $i++) {
        $index = $i + 1;
        $flag = 0;
        for (my $j = 0; $j <= $#pts; $j++) {
            if ($index == $pts[$j]) { $flag = $flag + 1; }
        }
        if ($flag != 0) { next; }
        print FILE $data[$i];
        $ptnum++;
    }
    close(FILE);
    print "New point number: $ptnum\n";
    return 0;
}

# remove_datapts_grad: remove data points from gradient file (natoms)
sub remove_datapts_grad {
    my ($flnm, $na, @pts) = @_;

    open(FILE, "<$flnm") or die "Could not open gradient file: $flnm!\n";
    my @data = grep(/\n/i, <FILE>);
    close(FILE);

    system("mv $flnm $flnm.old");

    my $index;
    my $rem;
    my $flag;
    my $ptnum = 0;
    open(FILE, ">$flnm");
    for (my $i = 0; $i <= $#data; $i++) {
        $index = $i / $na + 1;
        $flag = 0;
        for (my $j = 0; $j <= $#pts; $j++) {
            if ($index == $pts[$j]) { $flag = $flag + 1; }
        }
        if ($flag != 0) {
            $i = $i + $na - 1;
            next;
        }
        for (my $j = 0; $j < $na; $j++) {
            print FILE $data[$i + $j];
        }
        $i = $i + $na - 1;
        $ptnum++;
    }
    close(FILE);
    print "New point number: $ptnum\n";
    return 0;
}
