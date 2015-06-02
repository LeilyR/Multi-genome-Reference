#!/bin/bash
dir=/ebio/abt6_projects7/small_projects/mdubarry
fileReference=before_rr.fasta
fileRead=before_rr_read.fasta
fileAll=all.fasta

cd $dir/Documents/SampleProgram/bin/output
#rm * .d* dbReference.dbRead.las #rm all files
#Fasta prepare (format the ids)
#~/graph_git2/rungraph.pl fasta_prepare ./$fileReference $dir/Documents/pbsim/$fileReference
#~/graph_git2/rungraph.pl fasta_prepare ./$fileRead $dir/Documents/pbsim/$fileRead
cat $fileReference $fileRead > $fileAll

# Create Database
	#Reference
export PATH=$dir/software/DAZZ_DB:$PATH
fasta2DAM ./dbReference $fileReference
	#Read
fasta2DAM ./dbRead $fileRead

# Align read against reference
export PATH=$dir/software/DALIGNER:$PATH
HPCmapper -v $dir/Documents/SampleProgram/bin/output/dbReference $dir/Documents/SampleProgram/bin/output/dbRead | csh -v
#LAshow to look at the alignment
LAshow dbReference.dam dbRead.dam dbReference.dbRead.las -ca
# Convert .las to .sam using proovread
export PATH=$dir/software/proovread/bin:$PATH
export PERL5LIB=$dir/software/proovread/lib:$PERL5LIB
dazz2sam --las dbReference.dbRead.las --ref dbReference --qry dbRead > referenceRead.sam
#rm dbReference.dbRead.las #no need to keep this file, convert in sam format
