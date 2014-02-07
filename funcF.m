function F = funcF(dt, n, N_h0, e, xero)
    c = 0;     % integer
    N_h = N_h0 - dt*n;    % rad, new mean anomaly value
    while ( N_h >= 2*pi || N_h <= -2*pi )
        if ( N_h < 0 )
            N_h = N_h + 2*pi;   % bring down to -2pi < N < 0 rad
        else
            N_h = N_h - 2*pi;   % bring down to 0 < N < 2pi rad
        end
    end
    F = N_h / (e-1);      % rad, 1st approximation of new eccentric anomaly

    if ( N_h < (pi-xero) )
        N_h_try1 = e*sinh(F) - F;
        N_h_try2 = e*sinh(pi-F) - (pi-F);
        if ( abs(N_h_try2 - N_h) < abs(N_h_try1 - N_h) )
            F = pi-F;
        end
        % while value of N_h calculated from F is not-equal-to N_h...
        while ( abs( (e*sinh(F) - F) - N_h ) > 1E-8 ... % (tolerance = 1E-8)
                && abs( (e*(sinh(pi-F)) - (pi-F)) - N_h ) > 1E-8 ) 
            F_prev = F;     % save previous value of F
%            N_h_try = e*sinh(F) - F;
            F = F - (e*sinh(F) - F - N_h)/10; % rad, iterate new value of F
            c = c + 1;
            if ( abs( F - F_prev ) < 1E-10 )    % if it's going nowhere
                break;  % exit the while loop
            end
        end
    elseif ( N_h > (pi+xero) )
        N_h_try1 = e*sinh(F) - F;
        N_h_try2 = e*sinh(-pi-F) - (-pi-F);
        if ( abs(N_h_try2 - N_h) < abs(N_h_try1 - N_h) )
            F = -pi-F;
        end
        while ( abs( (e*(sinh(-pi-F)) + (pi+F)) - N_h ) > 1E-8 ...
                && abs( ( (e*sinh(F) - F) - N_h ) ) > 1E-8 )
            F_prev = F;     % save previous value of F
            
            F = F - (e*sinh(F) - F - N_h)/10;     % rad, iterate new value of F
            c = c + 1;
            if ( abs( F - F_prev ) < 1E-10 ) % if it's going nowhere
                break;  % exit the while loop
            end
        end
    else
        if ( N_h < xero && N_h > -xero )
            F = 0;
        else
            F = pi;
        end
    end    
    if ( abs( (e*sinh(F) - F) - N_h ) > xero )
        error('F did not converge');
    end
end
