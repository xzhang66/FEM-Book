function d=setup_ID_LM(d); 
include_flags;

count = 0; count1 = 0; 
for i = 1:neq
   if flags(i) == 2 % check if a node on essential boundary
      count = count + 1;
      ID(i) = count; % arrange essential B.C nodes first
      d(count)= e_bc(i); % store essential B.C in reordered form (d_bar)
   else
      count1 = count1 + 1;
      ID(i) = nd + count1;
   end 
end

for i = 1:nel 
   n = 1;
   for j = 1:nen
      blk = ndof*(IEN(j,i)-1); 
      for k = 1:ndof
         LM(n,i) = ID( blk + k );
         n = n + 1;
      end 
   end
end
