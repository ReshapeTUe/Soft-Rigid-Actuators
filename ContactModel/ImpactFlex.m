%% Detect impact between end-effector and ground
function [value,isterminal,direction] = ImpactFlex(~,q)
global l0


value = sin(q(2))*q(1)-sin(-0.25*pi)*l0; % Entry becomes zero is that contact point touches impact surface
isterminal = 1;     % Stop integration when impact occurs
direction = -1;   % Only detect impact if it occurs with negative velocity
end