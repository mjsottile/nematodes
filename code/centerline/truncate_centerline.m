function t = truncate_centerline(cline, startpt,row)

  [rs,cs] = size(cline);

  i = startpt;
  while (i <= rs && cline(i,1) > row)
      i = i+1;
  end
  
  if (i==rs)
      t = cline;
  else
      t = cline(1:i,:);
  end
  
  t = flipud(t);