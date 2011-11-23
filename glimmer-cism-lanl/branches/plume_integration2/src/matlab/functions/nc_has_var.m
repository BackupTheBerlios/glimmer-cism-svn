function hasVar = nc_has_var(f,varname)

  finfo = ncinfo(f);

  numVars = length(finfo.Variables);

  hasVar = false;
  for i=1:numVars
    thisVarName = finfo.Variables(i).Name;
  if (strcmp(thisVarName,varname))
    hasVar = true;
  end 

  end
end
