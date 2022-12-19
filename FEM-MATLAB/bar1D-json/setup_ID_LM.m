% setup ID and LM arrays
function d = setup_ID_LM(d);
include_flags;

count = 0; 
count1 = 0;

for i = 1:neq
   if flags(i) == 2 % check if essential boundary
      count = count + 1;
      ID(i) = count; % number first the nodes on essential boundary
      d(count)= e_bc(i); % store the reordered values of essential B.C
   else
      count1 = count1 + 1;
      ID(i) = nd + count1;
   end
end

for i = 1:nel
   for j = 1:nen
      LM(j,i)=ID(IEN(j,i)); % create the LM matrix
   end
end
