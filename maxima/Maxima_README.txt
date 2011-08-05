*******************************************************************************
README Notes for wxMaxima 
*******************************************************************************

1) At install, Maxima sets the enviroment variable 'maxima_userdir' to
   point at the curent user's '\maxima' directory - even if the directory
   doesn't exist.
   
   From the man pages:
   maxima_userdir:  Points to a directory for user customization files. 
                    Maxima's default search paths include maxima_userdir. 
				    Default value: $HOME/.maxima.
   
2) At install, Maxima sets the environment variable 'file_search_maxima' to
   include the directory specified in 'maxima_userdir' and to also include  
   '<maxima install root>\Maxima-<release #>\share\maxima\<release #>\share' 
   and all its subdirectories, recursed down to the lowest level.
  
   <maxima install root> - where the installer is pointed at install time
   (e.g. c:\users\conrad\applications\)
   
   Maxima-<release #> - concatentation of the string 'Maxima' and 
   the release value (e.g. Maxima-5.25.0)
   
   NOTE:  The file wild card is given by ###.{mac,mc}
   
3) At install, Maxima sets the environment variable 'file_search_lisp' to 
   include the directory specified in 'maxima_userdir' and to also include
   '<maxima install root>\Maxima-<release #>\share\maxima\<release #>\share'
   and all its subdirectories, recursed down to the lowest level, and also
   '<maxima install root>\Maxima-<release #>\src', pluse '$HOME/.lisp'.

   <maxima install root> - where the installer is pointed at install time
   (e.g. c:\users\conrad\applications\)
   
   Maxima-<release #> - concatentation of the string 'Maxima' and 
   the release value (e.g. Maxima-5.25.0)
   
   NOTE:  The file wild card is given by ###.{lisp,lsp}
   
4) The user can add additional directories to the environmental variables
   'file_search_maxima' and 'file_search_lisp' by 
   
   a)  creating a placing a file called 'maxima-init.mac' into 
       the 'maxima_userdir' location
   b)  with the lines 
           'file_search_maxima : append(["c:/cschiff/development/maxima/###.{mac,mc}"],file_search_maxima)$'
           'file_search_lisp   : append(["c:/cschiff/development/maxima/###.lisp"],    file_search_lisp)$'

5)  To get a package into wxMaxima
    example with the path c:\cschiff\Personal\projects\maxima\MyPackages and 
    the desired file 'critpts.mac' the command 

    load("C:/cschiff/Personal/projects/maxima/MyPackages/critpts.mac");

    will load 'critpts.mac' into Maxima memory.

	
6)  Line terminators:  ';' to end with an echo of the result to output and '$' to end with no output


    maxima_tempdir     : "c:/users/cschiff/maxima"$