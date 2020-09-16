function [fl,fid]=read_dlm_file(fname,dlm,nlines)
% READ_DLM_FILE  Read all or some of a delimited file into cell arrays-of-arrays.
% [FL,FID] = READ_DLM_FILE(FNAME,DLM,NLINES)
%    Reads a file (FNAME) and partitions to line using the delimeter (DLM).
%    The default DLM value is \t (9). The optional NLINES parameter
%    limits the number of lines read. 
%    
%    FL is a cell array of cell arrays of strings.  FL{n} is a cell array 
%    strings containing the nth line of the file.  
%    FID is the fid for the opened file, FNAME.
%
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

  if nargin==1
    dlm=9;
  end
  
  if ischar(fname)
    fid=fopen(fname,'r');
  else
    fid=fname;
  end
  
  ln=1;
  if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
  else
    do_close=0;
  end
  
  had_output=0;
  fl={};
  while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fl{ln}=dlmsep(tline,dlm);         
    ln=ln+1;
    if mod(ln,1000)==0
      verbose(['...' num2str(ln)],30);
      had_output=1;
    end
  end
  if do_close
    fclose(fid);
    fid=-1;
  end
  ln=ln-1;
  if had_output
    verbose([newline],30);
  end
end