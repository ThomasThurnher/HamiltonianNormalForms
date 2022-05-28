function [] = compare_flows(flow1,flow2,order)
'Comparing flows at normalisation order'
order

for i = 1:numel(flow1)
    flow1(i).coeffs
    flow2(i).coeffs

end
end