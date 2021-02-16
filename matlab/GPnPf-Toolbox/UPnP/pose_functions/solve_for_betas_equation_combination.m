function betas_and_f= solve_for_betas_equation_combination( beta_eq, betas_fij )

% -------------------------------------------------------------------------
% We obtain the value of the betas
% -------------------------------------------------------------------------
num_beta_eq= 12;
num_beta_unk= 3;
num_beta_unk_and_f= num_beta_unk + 1;
num_Bij= round(num_beta_eq/2);

% Due to the form the equations take we don't consider f for the
% linearization coordinates
linearization_coordinates= generate_linearization_coordinates(num_beta_eq/2, num_beta_unk);

% We will apply logarithms on both sides in order to calculate the
% absolute value
matrix_to_resolve= zeros(num_beta_unk_and_f);
independent_term= zeros(1,num_beta_unk_and_f);
for iterator=1:num_beta_unk_and_f
    % solution = [B1,B2, .. ,f]
    % If the equation is one of the Bfij
    if (beta_eq(iterator) > num_Bij)
        coordinates= linearization_coordinates(beta_eq(iterator)-num_Bij,:);
        % We introduce the f value for the vector
        matrix_to_resolve(iterator, num_beta_unk_and_f)= 2;
    else
        coordinates= linearization_coordinates(beta_eq(iterator),:);
    end

    if (coordinates(1)==coordinates(2))
        matrix_to_resolve(iterator, coordinates(1))= 2;
    else
        matrix_to_resolve(iterator, coordinates(1))= 1;
        matrix_to_resolve(iterator, coordinates(2))= 1;
    end

    independent_term(iterator)= log(abs(betas_fij(beta_eq(iterator))));
end

if (rank(matrix_to_resolve) == num_beta_unk_and_f)
    % betas_and_f= exp(inv(matrix_to_resolve)*independent_term');
    betas_and_f= exp(matrix_to_resolve\(independent_term'));

    for iterator=2:(num_beta_unk)
        betas_and_f(iterator)= betas_and_f(iterator)*sign(betas_fij(iterator));
    end
    betas_and_f(num_beta_unk_and_f)= abs(betas_and_f(num_beta_unk_and_f)); % f is always possitive
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create the vector to trasform the coordinates of the vector to the
% coordinates of the linearization
function linearization_coordinates= generate_linearization_coordinates(number_of_equations, number_of_unknowns)

coordinate_matrix= zeros(number_of_unknowns);
counter= 1;
for i=1:number_of_unknowns
    for j=i:number_of_unknowns
        coordinate_matrix(i,j)= counter;
        counter= counter+1;
    end
end

linearization_coordinates= zeros(number_of_equations,2);
for iterator=1:number_of_equations
    [row,column]= find(coordinate_matrix==iterator);
    linearization_coordinates(iterator,:)= [row, column];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%