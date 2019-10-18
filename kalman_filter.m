
function kalman_filter()

    load runway_accel_data.mat time_pts...
        std_dev_pos std_dev_spd std_dev_acc...
        s_meas v_meas a_meas

    % changes based on problem 
    n_states = 2;
    C = eye(n_states);
    dt = time_pts(2) - time_pts(1);
    R = [std_dev_pos^2 0; 0 std_dev_spd^2];
    G = std_dev_acc^2;
        
    % initialize 
    x_hat = [s_meas(:, 1); v_meas(:, 1)];
    P = R; 
    
    % plotting
    x_hat_rec = zeros(2, numel(time_pts));
    x_hat_rec(:, 1) = x_hat; 
    
    for n = 2:numel(time_pts)
        
        % A and B Matrics
        A = [1 dt; 0 1];
        B = [0; dt];
            
        Q = B*G*B';
        
        % predictive update
        u = a_meas(:, n-1);
        x_minus = A*x_hat + B*u; % we can replace this with rk4 step for att. kinematics
        P_minus = A*P*A' + Q;
        
        % Kalman gain
        K = P_minus*C' / (C*P_minus*C' + R);
        
        % measurement update
        z = [s_meas(:, n); v_meas(:, n)];
        x_hat = x_minus + K*(z - C*x_minus);
        P = (eye(n_states) - K*C)*P_minus; 
        
        % sum of diagonal elements 
        trace(P);

%         fprintf('X_hat = {%.4f}\n', x_hat)
        x_hat_rec(:, n) = x_hat;
    end
    
    figure(1);
    subplot(211);
    plot(time_pts, x_hat_rec(1, :), 'b-', time_pts, s_meas(:), 'r*');
    ylabel('Position (m)');
    subplot(212);
    plot(time_pts, x_hat_rec(2, :), 'b-', time_pts, v_meas(:), 'r*');
    ylabel('Speed (m/s)')
    xlabel('Time (s)')    

end






