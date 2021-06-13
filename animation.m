close all;

% Animation

% Simulation time
range = length(out.tout);

% Output data from simulnk model
theta1 = out.theta1_data.signals.values(1:range);
theta2 = out.theta2_data.signals.values(1:range);

% Position of cart
x = out.x_data.signals.values(1:range);
y = zeros(range);

% Position of first mass node
x_1 = x + l1.*sin(theta1);
y_1 = l2.*cos(theta1);

% Position of second mass node
x_2 = x +l1.*sin(theta1) + l2.*sin(theta2);
y_2 = l1.*cos(theta1) + l1.*cos(theta2);

% Plot double inverted pendulum
pend = plot(x_1(1),y_1(1),'ob',x_2(1),y_2(1),'og', [x(1) x_1(1)], [y(1) y_1(1)],'-b', [x_1(1) x_2(1)], [y_1(1) y_2(1)], '-g');
% Plot cart
cart_width = 0.3;
cart_height = 0.1;
cart = rectangle('Position',[x(1)-cart_width/2 y(1)-cart_height/2 cart_width cart_height]);
yline(-cart_height/2,'black');
xlim([-3 3])
ylim([-3 3])

% Timer
timer = 0;
text_label = text(0.45, 0.95, "Timer: " + num2str(timer), 'Units', 'normalized');

% Refresh coordinates
for k = 2:length(x)
    % Cart position
    delete(cart)
    cart = rectangle('Position',[x(k)-cart_width/2 y(k)-cart_height/2 cart_width cart_height]);
    % First node position
    pend(1).XData = x_1(k);
    pend(1).YData = y_1(k);
    % Second node position
    pend(2).XData = x_2(k);
    pend(2).YData = y_2(k);
    % First rode
    pend(3).XData = [x(k) x_1(k)];
    pend(3).YData = [y(k) y_1(k)];
    % Second rode
    pend(4).XData = [x_1(k) x_2(k)];
    pend(4).YData = [y_1(k) y_2(k)];
    % Timer
    delete(text_label)
    timer = timer + 0.01;
    text_label = text(0.45, 0.95, "Timer: " + num2str(timer), 'Units', 'normalized');
    % Refresh
    pause(0.01)
    drawnow
end