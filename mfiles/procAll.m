todoers = dir('ProcessCtd*.m');

for i=2:length(todoers);
  todoers(i).name
  %try
    eval(todoers(i).name(1:end-2));
 % catch
    
 % end
end;