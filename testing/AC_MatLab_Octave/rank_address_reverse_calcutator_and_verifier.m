clc
clear

for my_rank=1:128

  # Check the addressing of ranks via the subarray within arrayfun
  # In this case - user input
  # my_rank = 128
  sections = [  8,  4,  4  ];
  max_address_sections=sections(1)*sections(2)*sections(3);

  if (max_address_sections < my_rank)
    disp ("Rank is too large for sections")
  else
    
################################################################################
# Calculate the rank_section out of my_rank and sections[x,y,z]     
    zremainder = mod(my_rank, sections(1)*sections(2));
    if (zremainder == 0)
      rank_section = [ sections(1), sections(2), (my_rank - zremainder) / (sections(1)*sections(2)) ];   
    else
      rank_section(3) = (my_rank - zremainder) / (sections(1)*sections(2));
      
      yremainder = mod(zremainder, sections(1));
      if (yremainder == 0)
        rank_section = [ sections(1), (zremainder-yremainder)/sections(1), rank_section(3)+1 ];
      else
        rank_section = [ yremainder, (zremainder-yremainder)/sections(1)+1, rank_section(3)+1 ];
      endif
    endif
################################################################################

    rank_section;
    
    calculated_rank = (rank_section(3)-1)*sections(1)*sections(2) + (rank_section(2)-1)*sections(1) + rank_section(1);

    if (calculated_rank == my_rank)
      printf ("My_rank: %i was calculated successfully.\n", my_rank)
    else
      printf ("My_rank: %i was not calculated successfully.\n", my_rank)
    endif
  endif

endfor