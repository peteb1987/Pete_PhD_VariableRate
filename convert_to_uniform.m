 
function [data]=convert_to_uniform(t_in,data_in,t_grid_out)

   data(1:length(t_grid_out))=nan;
   t_inc=min(diff(t_grid_out));
   diff_t=find(t_in>0);
   t_index=0;
   data_index=1;
   for time=t_in(1):t_inc:t_in(end)
       t_index=t_index+1;
       while(t_in(data_index)<=time)
          data_index=data_index+1;
       end   
       data(t_index)=data_in(data_index-1);
      
   end 
   