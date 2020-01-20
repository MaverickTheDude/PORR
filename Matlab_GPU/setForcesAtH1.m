function Q1A = setForcesAtH1(q, input, ind, dq)
% Generates gravity force at first handle (H1) for all bodies. The order of
% bodies is from A to Z in a matrix form, i.e. size(Q1A) == [3,n]
% added: friction force for some cases. dq is last in order to faciliate
% backward compatibility
Nbodies = length(q); Q1A = zeros(3,Nbodies);
g = 9.80665; data = getS(q, input); 

for i = 1:Nbodies
    m = data(i).M1(1,1);
    Q1A(:,i) = data(i).S1c * [0; -m*g; 0];
    
    if isempty(input{1}(i).force)
       continue;   end % No applied forces
    
    for k = 1 : length(input{1}(i).force)
        switch input{1}(i).force(k).type
            case 'custom'
                force = input{1}(i).force;
                Q = [Rot(data(i).fi) * force(k).F; force(k).T];
                S1P = calcS(data(i).fi, force(k).s1P);
                Q1A(:,i) = Q1A(:,i) + S1P * Q;
            case 'control'
                force = input{1}(i).force;
                if ~isempty(force.F)
                    Q = [Rot(data(i).fi) * force(k).F(:,ind); force(k).T(ind)];
                    S1P = calcS(data(i).fi, force(k).s1P);
                    Q1A(:,i) = Q1A(:,i) + S1P * Q;
                end
                if ~isempty(force.T)
                    Q1A(3,i) = Q1A(3,i) + force.T(ind);
                end
            case 'dampingTrans'
                u = input{1}(i).H; % kierunek dzialania sily (+ skladowa obrotowa)
                c = input{1}(i).force(k).F; % 'oszukana' nazwa, zeby nie generowac pol
                % note: Czlonu 'j' nie bierzemy pod uwage, bo na
                % razie to zawsze jest ground body.
                Q1A(:,i) = Q1A(:,i) - data(i).S1c * u * c * dq(i);
            case 'dampingRot'
                H = input{1}(i).H;
                c = input{1}(i).force(k).F;
                Q1A(:,i)   = Q1A(:,i)   - H * c * dq(i);
                Q1A(:,i-1) = Q1A(:,i-1) + H * c * dq(i);
        end
    end
end

end

function S = calcS(fi, sLoc)
Om = [0, -1; 1, 0]; 
S = eye(3);
S(3,1:2) = (Om*Rot(fi)*sLoc).';
end