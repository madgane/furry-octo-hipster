
currentDIR = dir('.');
currentString = sprintf('%s',pwd);
addpath(currentString);

for iDirectory = 1:length(currentDIR)
   if currentDIR(iDirectory,1).isdir
       
       contWithDirectory = 1;
       if strcmp(currentDIR(iDirectory,1).name(1,1),'.')
           contWithDirectory = 0;
       end
       
       if contWithDirectory
           currentString = sprintf('%s\\%s',pwd,currentDIR(iDirectory,1).name);
           addpath(currentString);
       end

   end
end