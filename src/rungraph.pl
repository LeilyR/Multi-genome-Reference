#!/usr/bin/perl 
use strict;
use warnings;
use Cwd 'abs_path';

my @p = split /\//, abs_path($0);

my $dir = ''; 
for(my $i=0; $i<$#p; ++$i) {
	$dir.='/' unless $i==0;
	$dir.=$p[$i];
}

my $args ='';
for(my $i=0; $i<=$#ARGV; ++$i) {
	$args.=" $ARGV[$i]";
}

$ENV{'PATH'} = "$ENV{'PATH'}:$dir/../dazz_db:$dir/../daligner:$dir/..:$dir/../proovread/bin";

print "$ENV{'PATH'}\n";

my $cmd = "export LD_LIBRARY_PATH=$dir\n$dir/graph$args";
print "$cmd\n";
exec($cmd);


