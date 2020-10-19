#!/usr/local/bin/perl
if(($#ARGV != 1))
{
  print STDERR "Two arguments needed: macroname, njobs.\n";
  exit;
}
if ($ARGV[1] < 1) { print STDERR "The njobs argument has to be a number.\n"; exit; }
$njobs = int($ARGV[1]);
print STDERR "Submitting ".$njobs." jobs ...\n";
$maindir = `pwd`;
chomp($maindir);

$logdir = "/tmp/kincsesd/condor_logs";
$outdir  = `pwd`;
chomp($outdir);
print STDERR "Maindir = $maindir\n";
$macro = $ARGV[0];
print STDERR "Macro name: $macro \n";
$i=0;
while ($i<$njobs)
{
  $execfile="$maindir/Temp/run_${macro}_$i.sh";
  $jobfile="$maindir/Temp/job_${macro}_$i.job";
  if (-e "$maindir/Temp/${macro}_$i.log")
  {
    system("rm $maindir/Temp/${macro}_$i.*");
  }

  open(ROOT,">$execfile") || die "Could not open temp file\n";
  print ROOT "#!/bin/csh\n";
  print ROOT "cd $maindir\n";
  print ROOT "mkdir Temp/${macro}_$i/\n";
  print ROOT "cp exe/$macro.exe Temp/${macro}_$i/\n";
  print ROOT "cd Temp/${macro}_$i/\n";
  print ROOT "echo \"Environment: \"\n";
  print ROOT "echo \" \"\n";
  print ROOT "printenv\n";
  print ROOT "echo \" \"\n";
  print ROOT "./$macro.exe $i\n";
  print ROOT "rm ./$macro.exe\n";
  print ROOT "mv *.* $outdir/Data/.\n";
  print ROOT "cd $maindir\n";
  print ROOT "rm -rf  Temp/${macro}_$i/\n";
  close ROOT;
  chmod(0744,"$execfile");

  open(ROOT,">$jobfile") || die "Could not open temp file\n";
  print ROOT "Universe        = vanilla\n";
  print ROOT "\n";
  print ROOT "Executable      = $execfile\n";
  print ROOT "Priority        = +20\n";
  print ROOT "GetEnv          = True\n";
  print ROOT "\n";
  print ROOT "Initialdir      = $maindir/\n";
  print ROOT "Output          = $maindir/Logs/${macro}_${i}.out\n";
  print ROOT "Error           = $maindir/Logs/${macro}_${i}.err\n";
#  print ROOT "Log             = $logdir/Logs/${macro}_${i}.log\n";
  print ROOT "Log             = $maindir/Logs/${macro}_${i}.log\n";
  print ROOT "\n";
  print ROOT "+Job_Type       = \"Standard\"\n";
  print ROOT "\n";
  print ROOT "Queue\n";
  close ROOT;

  system("condor_submit $jobfile");
  print STDERR "Job #$i submitted.\n";

  sleep 0.1;
  $i=$i+1;
  if ($i>=$njobs) { last; }
}
close LIST;

