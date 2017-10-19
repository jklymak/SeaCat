st = argv();
name = st{1}
stname = st{2}
com = sprintf('%s= load(''%s'');', stname, name)
eval(com);
name_out = strrep(name, 'New', '');
com = sprintf('save(''-v7'', ''%s'', ''%s'');', name_out, stname)
eval(com);
