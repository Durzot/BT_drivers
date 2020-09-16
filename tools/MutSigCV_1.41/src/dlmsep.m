function res=dlmsep(s,d)
  if nargin==1
    d=9; % tab
  end

  pos=find(ismember(s,d));
  if ~isempty(pos)
    pos=[ 0 pos length(s)+1];
    for i=1:(length(pos)-1)
      res{i}=s((pos(i)+1):(pos(i+1)-1));
    end
  else
    res{1}=s;
  end
end