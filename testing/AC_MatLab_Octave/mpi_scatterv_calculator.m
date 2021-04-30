clc
clear

# Steering part of script below function!
# Amount of Processors provided must be dividable by 2 ### or better 4

function rtrn = mpi_scatter(size_mpi)
  
switch size_mpi
case 6 
  sections = [  3,  2,  1  ];
case 10 
  sections = [  5,  2,  1  ];
case 12  
  sections = [  3,  2,  2  ];
case 20  
  sections = [  5,  2,  2  ];
case 24  
  sections = [  4,  3,  2  ];
case 36  
  sections = [  4,  3,  3  ];
case 40  
  sections = [  5,  4,  2  ];
case 48  
  sections = [  4,  4,  3  ];
case 72  
  sections = [  6,  4,  3  ];
case 80  
  sections = [  5,  4,  4  ];
case 96 
  sections = [  6,  4,  4  ];
case 144
   sections = [  6,  6,  4 ];
case 216
   sections = [  6,  6,  6 ];

otherwise
 true_size  = size_mpi - mod(size_mpi, 2);
 sections             = [ 1, 1, 1 ];

 sw=0;
 while sections(1)*sections(2)*sections(3) <= true_size                                    
    sections(1) = sections(1) + 1;
    sw=1;
    if sections(1)*sections(2)*sections(3) <= true_size
      sections(2) = sections(2) + 1;
      sw=2;
    endif
    if sections(1)*sections(2)*sections(3) <= true_size
      sections(3) = sections(3) + 1;
      sw=3;
    endif
endwhile
sections(sw) = sections(sw) - 1;
end
domains=sections(1)*sections(2)*sections(3);
rtrn=domains;
endfunction

for i=2:2:256
  rtrn = mpi_scatter(i);
  printf("In: %i   Out: %i \n", i, rtrn)
endfor
