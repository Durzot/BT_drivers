function x=as_column(x)
  if size(x,2)>1
    x=x';
  end
end