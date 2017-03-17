function x = LHS(number,dimension)
    % Random permutation to create Latin Hypercube
    x = zeros(number,dimension);
    for i = 1:dimension
        x(:,i) = randperm(number)';
    end
    % Randomize within Latin Hypercube voxel.
    x = (x-0.5)/number;
    rn = rand(size(x));
    rn = (rn-0.5)/number;
    x = x + rn;
end