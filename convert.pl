my %hash;
my $file=$ARGV[0]; chomp $file; ### Gene List File
print "\nReading $file for converting Laevis genes to gene symbols\n";
open(F,"$file");
while(my $data = <F>)
{
 $data =~ s/^\s+|\s+$//g;
 chomp $data;
 if($data =~/Xelaev/)
 {}
 else
 {
  my @arr=split(/\./,$data);
  my $len=scalar(@arr);
  if(($len == 1)||($len == 2))
  {
   $hash{$arr[0]}=0;
  }
  if($len == 3)
  {
   my $var=$arr[0]."\.".$arr[1];
   $hash{$var}=0;
  }
 }
}
close F;

my @arrf=split(/\.txt/,$file);
my $name=$arrf[0]."\_"."for"."\_"."ortho";
open(OUT,">$name.txt");
print "\n$file is processed and writing final list to $name\.txt\n";
foreach my $m(keys %hash)
{
 if($m eq "")
 {}
 else
 {
  print "$m\n";
  print OUT "$m\n";
 }
}
close OUT;
print "\nAll Done and Dusted\n";
print "\nPlease Contact Praneet Chaturvedi for issues and bugs\n";