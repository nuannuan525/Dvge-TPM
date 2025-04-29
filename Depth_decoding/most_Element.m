function [mostFrequentElement] = most_Element(idx_t1)

nonzero_t1 = idx_t1(idx_t1 ~= 0);
if isempty(nonzero_t1)
        mostFrequentElement=nan;
else
    uniqueElements = unique(nonzero_t1);
    counters = zeros(size(uniqueElements));

    % count the number of each elements
    for i = 1:length(uniqueElements)
        counters(i) = sum(nonzero_t1 == uniqueElements(i));
    end

    % Find the most frequently occurring element
    [~, index] = max(counters);
    mostFrequentElement = uniqueElements(index);

end

