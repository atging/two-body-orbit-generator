function E = funcE(dt, tau, M0, e, xero, oType) 
c = 0;     % integer
swtch = 1;
if ( oType <= 1 )   % Circular or Elliptical
    M = 2*pi*dt/tau + M0*pi/180;    % rad, new mean anomaly value
    while ( M >= 2*pi )
        M = M - 2*pi;   % bring down to 0 < M < 2pi rad
    end
    E = M / (1-e);      % rad, 1st approximation of new eccentric anomaly
            
    if ( M < (pi-xero) )
        M_try1 = E - e*sin(E);
        M_try2 = (pi-E) - e*sin(pi-E);
        if ( abs(M_try2 - M) < abs(M_try1 - M) )
            E = pi-E;
        end
        % while value of M calculated from E is not-equal-to M...
        while ( abs( (E - e*sin(E)) - M ) > 1E-8 ... % (tolerance = 1E-8)
                && abs( ((pi-E) - e*(sin(pi-E))) - M ) > 1E-8 ) 
            E_prev = E;     % save previous value of E
            if ( c > 15 && abs( (E - e*sin(E)) - M ) > 0.2 )
                % if haven't found E giving calculated M within 0.1 rad...
                swtch = -swtch;
                c = -10;
            end
            if ( swtch == 1 )
                E = M + e*sin(E);     % rad, iterate new value of E
            else    % must be alternate value of sin(E)
                E = M + e*sin(pi-E);    % rad
            end
            c = c + 1;
            if ( abs( E - E_prev ) < 1E-10 )    % if it's going nowhere
                break;  % exit the while loop
            end
        end
    elseif ( M > (pi+xero) )
        M_try1 = E - e*sin(E);
        M_try2 = (-pi-E) - e*sin(-pi-E);
        if ( abs(M_try2 - M) < abs(M_try1 - M) )
            E = -pi-E;
        end
        while ( abs( (-(pi+E) - e*(sin(-pi-E))) - M ) > 1E-8 ...
                && abs( ( (E - e*sin(E)) - M ) ) > 1E-8 )
            E_prev = E;     % save previous value of E
            if ( c > 15 && abs( (E - e*sin(E)) - M ) > 0.1 )
                % if haven't found E giving calculated M within 0.1 rad...
                swtch = -swtch;
                c = -10;
            end
            if ( swtch == 1 )
                E = M + e*sin(E);     % rad, iterate new value of E
            else    % must be alternate value of sin(E)
                E = M + e*sin(-pi-E);    % rad
            end
            c = c + 1;
            if ( abs( E - E_prev ) < 1E-10 )    % if it's going nowhere
                break;  % exit the while loop
            end
        end
    else
        if ( M < xero && M > -xero )
            E = 0;
        else
            E = pi;
        end
    end    
end