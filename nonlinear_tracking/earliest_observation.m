function [ ind ] = earliest_observation( list, time )
%EARLIEST_OBSERVATION Find the index of the earliest observation after a
%particular changepoint

ind = find(list==min(list(list>time)),1);

end

