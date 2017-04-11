function h = get_circular_filter(rad)
    h = zeros(2*rad+1,2*rad+1);
    [X,Y] = meshgrid(-rad:rad,-rad:rad);
    mask = abs(X).^2 +abs(Y).^2 <= rad^2;
    h(mask) = 1;
end