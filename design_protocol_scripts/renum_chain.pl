#!/usr/bin/perl
my $n = 0;
my $old_count = "";
my $old_chain = "";
if($ARGV[0]){$n = $ARGV[0];}
# ATOM   9953 3HG2 THR A 202     -13.364  17.786  10.102  1.00  0.00           H
while($line = <STDIN>)
{
	chomp($line);
	if(substr($line,0,6) eq "HETATM" and substr($line,17,3) eq "MSE")
	{
		$line =~ s/^HETATM/ATOM  /;
		substr($line,17,3,"MET");
		substr($line,12,2,"SD") if substr($line,12,2) eq "SE";
	}
	if(substr($line,0,4) eq "ATOM")
	{
		$new_count = substr($line,22,4);
		$new_chain = substr($line,21,1);
		if($old_chain eq "" or $new_chain ne $old_chain){$n = 0; $old_count = ""; $old_chain = $new_chain;}
		if($old_count eq "" or $new_count != $old_count){$n++;$old_count = $new_count; }
		
		print substr($line,0,22) . sprintf("%4s",$n) . substr($line,26)."\n";
	}
}
