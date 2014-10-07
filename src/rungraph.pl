#!/usr/bin/perl 
use strict;
use warnings;

my @p = split /\//, $0;

my $dir = ''; 
for(my $i=0; $i<$#p; ++$i) {
	$dir.='/' unless $i==0;
	$dir.=$p[$i];
}

my $args ='';
for(my $i=0; $i<=$#ARGV; ++$i) {
	$args.=" $ARGV[$i]";
}

my $cmd = "export LD_LIBRARY_PATH=$dir\n$dir/graph$args";

print "$cmd\n";

exec($cmd);



