file = 'solfield.dat';

out = fopen('finverse.dat','w');

M = load(file);

I = pinv(M);

for i=1:length(I);
   for j=1:length(I(i,:));
      fprintf(out,'%f\t',I(i,j));
   end
   fprintf(out,'\n');
end

fclose(out);
