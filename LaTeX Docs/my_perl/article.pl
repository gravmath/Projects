#!/usr/bin/perl
#
#  article.pl
#
use Tk;

sub text_translate()
{
  
  use translate;

  $strc_filename = $en1_fr1->get();
  $out_filename  = $en1_fr2->get();  
  $title         = $en1_fr4->get();
  $author        = $en1_fr5->get();
  
 open(STRC, $strc_filename);
 open(OUT,  ">$out_filename");

 @strc = ();
 
 while( <STRC> )
 {
   push(@strc,$_);
 }

 print OUT "\\documentclass\[12pt\]\{article\}\n";
 print OUT "\\usepackage\{latexsym\}\n";
 print OUT "\\usepackage\{epic,eepic,graphicx,url,amsmath,amssymb,latexsym,alltt\}\n";
 print OUT "\\begin\{document\}\n";
 print OUT "\\title\{$title\}\n";
 print OUT "\\author{$author}\n";
 print OUT "\\maketitle\n\n";

 foreach $strc (@strc)
 {
   if( $strc =~ /dat$/ )
   {
     open(IN, $strc) or die "Can't open $strc";
     while ( <IN> )
     {
       @words = split(" ",$_);
       foreach $word (@words)
       {
         translate::translate(\$word);
         push(@twords,$word);
       }
       $line = join(" ",@twords);
       @twords = ();
       print OUT $line, "\n";
     }
     close(IN);
   }
   else
   {
     print OUT $strc, "\n";
   }
 }
 print OUT "\\end\{document\}";
 close(STRC);
 close(OUT);
}


sub  TexIt()
{
   $out_filename  = $en1_fr2->get();  
   system("latex $out_filename");

}
#***************************************************
#Define the main window
#***************************************************
$mw = MainWindow->new();
$mw->configure(-title      => 'LaTeX Interpretter',
              );
$mw->geometry('+20+20');  #position the upper left corner            

#***************************************************
#Define the first frame and its widgets
#***************************************************
$fr1  = $mw->Frame(-relief      => 'groove',
                      -borderwidth => 3,
                     )->pack(-side => 'top',
                             -fill => 'x',
                            );
                            
$lb1_fr1 = $fr1->Label(-text => 'Enter Section Filename     ',
                         )->pack(-side => 'left',
                                );
                         
$en1_fr1 = $fr1->Entry(-width => 50,
                         )->pack(-side => 'left',
                                 -pady => 10,
                                );
                     
#***************************************************
#Define the second frame and its widgets
#***************************************************
$fr2   = $mw->Frame(-relief       => 'groove',
                       -borderwidth  => 3,
                      )->pack(-side  => 'top',
                              -fill  => 'x',
                             );
                             
$lb1_fr2 = $fr2->Label(-text => 'Enter the LaTeX Filename',
                         )->pack(-side => 'left',
                                );

$en1_fr2 = $fr2->Entry(-width => 50,
                      )->pack(-side => 'bottom',
                              -fill => 'x',
                              -pady => 10,
                             );

#***************************************************
#Define the fourth frame and its widgets
#***************************************************
$fr4   = $mw->Frame(-relief       => 'groove',
                       -borderwidth  => 3,
                      )->pack(-side  => 'top',
                              -fill  => 'x',
                             );
                                                            
$lb1_fr4 = $fr4->Label(-text => 'Title                               ',
                         )->pack(-side => 'left',
                                );

$en1_fr4 = $fr4->Entry(-width => 50,
                      )->pack(-side => 'right',
                              -fill => 'x',
                              -pady => 10,
                             );
                             
#***************************************************
#Define the fifth frame and its widgets
#***************************************************
$fr5   = $mw->Frame(-relief       => 'groove',
                       -borderwidth  => 3,
                      )->pack(-side  => 'top',
                              -fill  => 'x',
                             );                             
$lb1_fr5 = $fr5->Label(-text => 'Author                               ',
                         )->pack(-side => 'left',
                                );

$en1_fr5 = $fr5->Entry(-width => 50,
                      )->pack(-side => 'right',
                              -fill => 'x',
                              -pady => 10,
                             );
                             
                             
                             

#***************************************************
#Define the third frame and its widgets
#***************************************************
$fr3   = $mw->Frame(-relief       => 'groove',
                       -borderwidth  => 3,
                      )->pack(-side  => 'top',
                              -fill  => 'x',
                             );

$exit_b_fr3 = $fr3->Button(-text    => "Exit",
                             -command => sub{ exit},
                             -relief  => groove,
                            )->pack(-side => 'right',
                                   );

$trans_b_fr3  = $fr3->Button(-text    => "Translate the File",
                             -command  => \&text_translate,
                             -relief   => groove,
                            )->pack(-side => 'left',
                                   );
                                   
$LaTex_b_fr3  = $fr3->Button(-text    => "LaTeX the File",
                             -command => \&TexIt,
                             -relief  => groove,
                            )->pack(-side => 'left',
                                    -padx => 100,
                                   );
                          
                                   

MainLoop;