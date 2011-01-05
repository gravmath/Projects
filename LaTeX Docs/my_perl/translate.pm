package translate;

#       "#"      => " ",
#       "#pro",   => "\\vspace\{5mm\}\n\\hrule\n\\begin\{alltt\}",
#       "#end",   => "\\end\{alltt\}\n\\hrule\n\\vspace\{5mm\}",
#
#       "#"      =>  "%",
#       "#pro",   => "%",
#       "#end",   => "%",

%h = ( 
       "#"      =>  " ",
       "#pro",   => "\\clearpage\n\\vspace\{5mm\}\n\\hrule\n\\begin\{alltt\}",
       "#end",   => "\\end\{alltt\}\n\\hrule\n\\vspace\{5mm\}",       
       "#a"      => "\\alpha",
       "#A"      => "\\Alpha",
       "#b"      => "\\beta",
       "#B"      => "\\Beta",
       "#c"      => "\\chi",
       "#C"      => "\\Chi",
       "#d"      => "\\delta",
       "#D"      => "\\Delta",
       "#e"      => "\\epsilon",
       "#E"      => "\\Epsilon",
       "#f"      => "\\phi",
       "#F"      => "\\Phi",
       "#g"      => "\\gamma",
       "#G"      => "\\Gamma",
       "#h"      => "\\eta",
       "#H"      => "\\Eta",
       "#i"      => "\\iota",
       "#I"      => "\\Iota",
       "#j"      => "\\epsilon",
       "#J"      => "\\Epsilon",
       "#k"      => "\\kappa",
       "#K"      => "\\Kappa",
       "#l"      => "\\lambda",
       "#L"      => "\\Lambda",
       "#m"      => "\\mu",
       "#M"      => "\\Mu",
       "#n"      => "\\nu",
       "#N"      => "\\Nu",
       "#o"      => "o",
       "#O"      => "\\Omicron",
       "#p"      => "\\pi",
       "#P"      => "\\Pi",
       "#q"      => "\\theta",
       "#Q"      => "\\Theta",
       "#r"      => "\\rho",
       "#R"      => "\\Rho",
       "#s"      => "\\sigma",
       "#S"      => "\\Sigma",
       "#t"      => "\\tau",
       "#T"      => "\\Tau",
       "#u"      => "\\upsilon",
       "#U"      => "\\Upsilon",
       "#v"      => "\\varpi",
       "#V"      => "\\Epsilon",
       "#w"      => "\\omega",
       "#W"      => "\\Omega",
       "#x"      => "\\xi",
       "#X"      => "\\Xi",
       "#y"      => "\\psi",
       "#Y"      => "\\Psi",
       "#z"      => "\\zeta",
       "#Z"      => "\\Zeta",
       "#lp"     => "\\left\(",
       "#rp"     => "\\right\)",
       "#lb"     => "\\left\[",
       "#rb"     => "\\right\]",
       "#lc"     => "\\left\\\{",
       "#rc"     => "\\right\\\}",
       "#ld"     => "\\left.",
       "#rd"     => "\\right.",
       "#ubo"    => "\\underbrace\{",
       "#bc"     => "\}",
       "#be"     => "\\begin\{equation\}\\label\{c0\}",
       "#ee"     => "\\end\{equation\}",
       "#bea"    => "\\begin\{eqnarray\}\\label\{c0\}",
       "#eea"    => "\\end\{eqnarray\}",
       "#ae"     => "& = &",
       "#as"     => "&   &",
       "#nn"     => "\\nonumber",
       "#ep"     => "\\quad \.",
       "#ec"     => "\\quad \,",
       "#ref"    => "Eq. \(\\ref\{c0\}\)",
       "#cite"   => "\[\\cite\{c0\}\]",
       "#box"    => "\\boxed\{",
       "#def"    => "\\equiv",
       "#LHS"    => "left-hand-side",
       "#RHS"    => "right-hand-side",
       "#bfig"   =>  "\\begin\{figure\} \n \\centerline\{ \n\t\\includegraphics\[height= c0 in,width= c1 in,keepaspectratio\]\{c2\}\}",
       "#efig"   =>  "\\label\{c0\} \n \\end\{figure\}",
       "#mc"     => "\\mathcal\{ c0 \}",
       "#nin"    => "\\noindent",
       "#in"     => "\\indent",
       "#vec"    => "\{\\vec c0\}",
       "#form"   => "\{\\tilde c0\}",
       "#bv"     => "\{\\vec c0\}_\{c1\}",
       "#bf"     => "\{\\tilde c0\}^\{c1\}",
       "#uv"     => "\{\\hat c0\}_\{c1\}",
       "#hat"    => "\{\\hat c0\}",
       "#til"    => "\{\\tilde c0\}",
       "#ip"     => "\\left\( c0 , c1\\right\)",
       "#dp"     => "c0 \\cdot c1",
       "#cp"     => "\\langle c0 , c1 \\rangle",
       "#Dp"     => "\\langle c0 | c1 \\rangle",
       "#SM2"    =>  "\\left\( \\begin\{array\}\{cc\}   \n c0 \& c1             \\\\ \n c2 \& c3                   \n                     \\end\{array\} \\right\)",
       "#SM3"    =>  "\\left\( \\begin\{array\}\{ccc\}  \n c0 \& c1 \& c2       \\\\ \n c3 \& c4 \& c5       \\\\  \n  c6 \& c7 \& c8  \n \\end\{array\} \\right\)",
       "#SM4"    =>  "\\left\( \\begin\{array\}\{cccc\} \n c0 \& c1 \& c2 \& c3 \\\\ \n c4 \& c5 \& c6 \& c7 \\\\\ \n  c8 \& c9 \& c10 \& c11 \\\\ c12 \& c13 \& c14 \& c15 \n \\end\{array\} \\right\)",
       "#ca2"    =>  "\\left\( \\begin\{array\}\{c\} \n c0 \\\\ c1 \\end\{array\} \\right\)",
       "#ca3"    =>  "\\left\( \\begin\{array\}\{c\} \n c0 \\\\ c1 \\\\ c2 \\end\{array\} \\right\)",
       "#ca4"    =>  "\\left\( \\begin\{array\}\{c\} \n c0 \\\\ c1 \\\\ c2 \\\\ c3 \\end\{array\} \\right\)",
       "#ca6"    =>  "\\left\( \\begin\{array\}\{c\} \n c0 \\\\ c1 \\\\ c2 \\\\ c3 \\\\ c4 \\\\ c5 \\end\{array\} \\right\)",       
       "#ra2"    =>  "\\left\( \\begin\{array\}\{cc\}     \n c0 \& c1                         \\end\{array\} \\right\)",       
       "#ra3"    =>  "\\left\( \\begin\{array\}\{ccc\}    \n c0 \& c1 \& c2                   \\end\{array\} \\right\)",              
       "#ra4"    =>  "\\left\( \\begin\{array\}\{cccc\  } \n c0 \& c1 \& c2 \& c3             \\end\{array\} \\right\)",              
       "#ra6"    =>  "\\left\( \\begin\{array\}\{cccccc\} \n c0 \& c1 \& c2 \& c3 \& c4 \& c5 \\end\{array\} \\right\)",                     
       "#dby"    => "\\frac\{ d c0\}\{d c1\}",
       "#cby"    => "{ c0 }_{ , c1 }",
       "#pby"    => "\\frac\{\\partial c0\}\{\\partial c1\}",
       "#vby"    => "\\frac\{\\delta c0\}\{\\delta c1\}",
       "#par"    => "\\partial_\{c0\}",
       "#var"    => "\\delta",
       "#fvar"   => "\{ \\delta c0 \} \|_\{\\delta c1 \}",
       "#int3"   => "\\int c0 r^2 \\sin\(\\theta\) dr d\\theta d\\phi",
       "#dfunc"  => "\\delta \\left\( c0 \\right\)",
       "#kd"     => "\\delta_\{c0 c1\}",
       "#Kd"     => "\{\\delta ^ c0\}_\{c1\}",
       "#Kdt"    => "\{\\delta _ c0\}^\{c1\}",
       "#KD"     => "\\delta^\{c0 c1\}",
       "#cint"   => "\\int d^3x dt",
       "#aint"   => "\\int d^3a dt",
       "#lint"   => "\\int d \\lambda",
       "#Tu"     => "\{ c0 \} ^ \{ c1 \}",
       "#Td"     => "\{ c0 \} _ \{ c1 \}",       
       "#Tud"    => "\{ \{ c0 } ^ \{ c1 } \}_\{ c2 \}",
       "#Tdu"    => "\{ \{ c0 } _ \{ c1 } \}^\{ c2 \}",
       "#symo"   => "\(c0",
       "#symc"   => "c0\)",
       "#asymo"  => "\[c0",
       "#asymc"  => "c0\]",
       "#3p1"    => "3+1",
       "#1_16p"  => "\\frac\{c0\}\{16 \\pi\}",
       "#1_2"    => "\\frac{c0}{2}",
       "#nu"     => "\{n ^ c0\}",
       "#nd"     => "\{n _ c0\}",
       "#proju"  => "\{\\perp \} ^ \{ c0 c1 \}",
       "#proj"   => "\{\\perp ^ c0}_\{ c1 \}",
       "#projt"  => "\{\\perp _ c0}^\{ c1 \}",
       "#projd"  => "\{\\perp \} _ \{ c0 c1 \}",
       "#lapse"  => "\\alpha",
       "#shift"  => "\{\\vec \\beta\}",
       "#shiftu" => "\{\\beta\}^\{c0\}",
       "#shiftd" => "\{\\beta\}_\{c0\}",
       "#3g"     => "\{\\gamma\}_\{c0 c1\}",
       "#3ig"    => "\{\\gamma\}^\{c0 c1\}",
       "#3gm"    => "\{\\gamma ^c0 \}_\{ c1 \}",
       "#r3g"    => "\\sqrt{\\gamma\}",
       "#det3g"  => "\\gamma",
       "#4g"     => "\{g\}_\{c0 c1\}",
       "#4ig"    => "\{g\}^\{c0 c1\}",
       "#Tran"   => "\{\\Lambda\}^\{c0\}_\{c1\}",
       "#Trant"  => "\{\\Lambda\}_\{c0\}^\{c1\}",
       "#Cu"     => "\{\\nabla\}^\{c0\}",
       "#Cd"     => "\{\\nabla\}_\{c0\}",
       "#Du"     => "D^\{c0\}",
       "#Dd"     => "D_\{c0\}",
       "#lie"    => "\\pounds_{c0} c1",
       "#Chr"    => "\{ \\Gamma ^\{c0\}\}_\{c1 c2\}",
       "#3Chr"   => "\{ ^\{\(3\)\} \\Gamma ^\{c0\}\}_\{c1 c2\}",
       "#Riem"   => "\{R\}^\{c0\}_\{c1 c2 c3\}",
       "#3R"     => "R",
       "#3Ru"    => "R^\{c0 c1\}",
       "#3Rd"    => "R_\{c0 c1\}",
       "#3Rm"    => "\{R^c0\}_\{c1\}",
       "#R0"     => "R^0",
       "#Ri"     => "R^\{c0\}",
       "#Qa"     => "\\alpha R^0",
       "#Qb"     => "\\beta_{c0} R^{c1}",
       "#pi"     => "\{\\pi\}^\{c0 c1\}",
       "#pim"    => "\{\\pi^c0\}_\{c1\}",
       "#pid"    => "\{\\pi\}_\{c0 c1\}",
       "#TrP"    => "Tr\( \\pi \)",
       "#TrP2"   => "Tr\( \\pi ^2 \)",       
       "#3K"     => "\{K\}^\{c0 c1\}",
       "#3Km"    => "\{K^c0\}_\{c1\}",
       "#3Kd"    => "\{K\}_\{c0 c1\}",
       "#TrK"    => "Tr\( K \)",
       "#TrK2"   => "Tr\( K ^2 \)",
       "#zpos"   => "{z ^ c0}",
       "#zdot"   => "\{\\dot z\}^\{c0\}",
       "#uvel"   => "\{u\}_\{c0\}",
       "#uvelu"  => "\{u\}^\{c0\}",
       "#unorm"  => "|| u ||",
       "#muu"    => "\\sqrt\{1 + ||u||^2\}",
       "#ep"     => "\{\\epsilon \}^\{ c0 \}",
       "#enorm"  => "|| \\epsilon ||",
       "#mue"    => "\\sqrt\{1 - || \\epsilon || ^2 / \\alpha ^2 \}",
       "#Ham"    => "\\mathcal\{H\}",
       "#Jac"    => "\{J\}^\{c0\}_\{c1\}",
       "#Jact"   => "\{J\}_\{c0\}^\{c1\}",
       "#detJ"   => "|J|",       
       "#rho0"   => "\{ \\tilde \\rho \}_\{0\}",
       "#hh"     => "\( 1 + e \)",
       "#Schw_f" => "\\left\( 1 - \\frac\{2M\}\{r\} \\right\)",
       "#W2"     => "\\frac\{105\}\{32 \\pi h^3 \} \\left\(1 - \\frac\{|\\vec x - \\vec z|^2\}\{h^2\} \\right\)^2",
       "#W3"     => "\\frac\{315\}\{64 \\pi h^3 \} \\left\(1 - \\frac\{|\\vec x - \\vec z|^2\}\{h^2\} \\right\)^3",
       "#W4"     => "\\frac\{3465\}\{512 \\pi h^3 \} \\left\(1 - \\frac\{|\\vec x - \\vec z|^2\}\{h^2\} \\right\)^4",
       "#Wa"     => "W \\left\( |\\vec x - \\vec z| ; h \\right\)");



#*********************************************************
sub translate
{
  my ($input) = @_;
  my @a = ();
  my @b = ();
  my @c = ();
  my $c;

  #step 1 - if there is evidence of a command isolate the 
  #         command and arguments - recurse through the 
  #         arguments in the event there is a nested command 
  if( $$input =~ /#.{1,}\(/ )
  {
    @a = find_command($$input);
    @b = find_args($a[1]);
    foreach $b (@b)
    {
      if ( $b =~ /#.{1,}\(/ )
      {
        translate(\$b);
      }
      else
      {
        if ( $b =~ /#/ )
        {
          $b =~ s/$b/$h{$b}/;
        }  
      }
      push(@c,$b);    
    }
    $command = $h{$a[0]};
    $i = 0;
    foreach $c (@c)
    {
      $arg = "c$i";
      $command =~ s/$arg/$c[$i]/;
      $i++;
    }
    $$input = $command;
  }
  else
  {
    if ( $$input =~ /#/ )
    {
      $$input =~ s/$$input/$h{$$input}/;
    }
  }
}

#*********************************************************
sub find_command
{
  my ($input) = @_;
  my @a       = split("",$input);
  my $base    = 0;
  my $size    = @a;
  my $counter = 0;
  my @ret     = ();
  my $end     = $size - 1;
  
  #count the characters until the first open parenthesis is encountered 
  foreach $a (@a)
  {
    #first check to see if the current character matches an
    #opened parenthesis
    if ( $a =~  /\(/ )
    {
      last;
    }  
    $counter += 1;
  }
    
  push(@ret,substr($input,$base,$counter));
  push(@ret,substr($input,$counter+1,$end-$counter-1));
  
  return @ret; 
}


#*********************************************************
sub find_args
{
  my ($input) = @_;
  my @a       = split("",$input);
  my $op      = 0;
  my $cp      = 0;
  my $bal     = 0;
  my $counter = 0;
  my @loc     = ();
  my @args    = ();
  my $num_com = 0;
  my $cur_loc = 0;
  my $nxt_loc = 0;
  my $size    = @a;

  #count the number of 'outside' commas and note their locations  
  foreach $a (@a)
  {
    #first check to see if the current character matches an
    #open or closed parenthesis or a comma
    if ( $a =~  /\(/ )
    {
      $op += 1;
    }  
    if ( $a =~ /\)/ )
    {
      $cp += 1;
    }
    if ( $a =~ /,/  && ( $op - $cp == 0 ) )
    {
      $num_com += 1; 
      push(@loc,$counter);
    }
    $counter += 1;
  }
  if( $num_com > 0)
  {
    #print "The array \@loc is:", @loc,"\n";
  
    #first take care of the initial boundary
    $cur_loc = 0;
    $nxt_loc = $loc[0];
    #print "first: $cur_loc, $nxt_loc ",substr($input,$cur_loc,$nxt_loc-$cur_loc),"\n";
    push(@args,substr($input,$cur_loc,$nxt_loc-$cur_loc));

    #next do the interior
    for ( $i = 1; $i < $num_com; $i++)
    {
      $cur_loc = $nxt_loc +1;
      $nxt_loc = $loc[$i];
      #print "$cur_loc, $nxt_loc ",substr($input,$cur_loc,$nxt_loc-$cur_loc),"\n";
      push(@args,substr($input,$cur_loc,$nxt_loc-$cur_loc));
    }
  
    #finish with the last boundary
    $cur_loc = $nxt_loc + 1;
    $nxt_loc = $size;
    #print "last: $cur_loc, $nxt_loc ",substr($input,$cur_loc,$nxt_loc-$cur_loc),"\n";
    push(@args,substr($input,$cur_loc,$nxt_loc));
  }
  else
  {
     push(@args,$input);
  }
  return @args;
}

1;