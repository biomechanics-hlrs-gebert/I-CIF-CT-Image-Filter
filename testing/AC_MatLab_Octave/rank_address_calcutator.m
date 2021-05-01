clc
clear

# Check the addressing of ranks via the subarray within arrayfun
# In this case - user input
sections = [  8,  4,  4  ];

domains=sections(1)*sections(2)*sections(3)

testarray = zeros(1,domains)

disp("")

counter=0;
for ii=1:sections(1)
  for jj=1:sections(2)
    for kk=1:sections(3)
        address = (kk-1)*sections(1)*sections(2) + (jj-1)*sections(1) + ii
        testarray (address) = 1;
        counter = counter +1;
    endfor
  endfor
endfor

disp("")

counter

disp("")

testarray