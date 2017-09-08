%%% Code to do Cobwebbing for random regular graphs in the thermodynamic limti


% Parameters
lambda = 3.4379; % Tested lambda
A_0 = rand(); % Specify the initial value of A
K = 4; % Degree
J = 1;

pts = 60; % number of points of the linspace where the function is going to be plotted
steps = 20; % number of steps that are going to be taken
lim1 = 5; % limits of the linspace
lim2 = .2;
%---------------------------------------------------

% Set up two curves to create a plot for cobwebbing:
As1 = linspace(-lim1, -lim2, pts);
As2 = linspace(lim2, lim1, pts); % Two linspaces are done bechause it is a hyperbola

fAs1 = lambda-(K-1)*J^2./As1;
fAs2 = lambda-(K-1)*J^2./As2;  % Function that wants to be tested

%---------------------------------------------------

figure % This creates a new graphics window
plot(As1,fAs1,'k-', 'DisplayName', '\lambda=3')
hold on

plot(As2, fAs2, 'k-');% Make the plot
%hold on
xlabel('A^*') % Add some labels
ylabel('f(A^*)');
legend('\lambda=3.4');
grid on
%string_for_title = ['Recurrence relation for \lambda = ',num2str(lambda)];
%title(string_for_title)
% Now add a line with slope = 1

plot([-lim1 lim1],[-lim1 lim1],'k:', 'LineWidth', 1.2)
plot([0,0], [-5, 7], 'k:', 'LineWidth', 1.2)
plot([-5,5], [0, 0], 'k:', 'LineWidth', 1.2)

%---------------------------------------------------

% This is how you write the contents of a graphics window to a file (type
% "help print" for additional details:
print -depsc cob0.eps
% If you want a different format, you could comment out the above and use
% one of the below instead:
% print -dpng cob1.png
% print -djpeg100 cob1.jpg

%---------------------------------------------------

% Put in the first cobweb line. This is a line that goes from the point
% P_initial on the bottom of the plot, up to the corresponding point on
% the recurrence relation:
x_coords = [A_0 A_0];
A_1 = lambda-(K-1)*J^2/A_0;
y_coords = [0 A_1];
plot(x_coords,y_coords,'k--')
% This is how you place text labels on the plot:
%text(mean(x_coords) - 0.035,mean(y_coords),’#1’)
% Save the current graphics window to a file
print -depsc cob1.eps

%---------------------------------------------------

for i = 1:steps

    % Put in the second cobweb line. This is a line that goes horizontally
    % from the previous point to the line of slope 1.
    x_coords = [A_0 A_1];
    y_coords = [A_1 A_1];
    plot(x_coords,y_coords,'k--')
    % Place text label on the plot:
    %text(mean(x_coords),mean(y_coords) + 0.035,’#2’)
    % Save the current graphics window to a file
    print -depsc cob2.eps

    A_0 = A_1;

    % Put in the first cobweb line. This is a line that goes from the point
    % P_initial on the bottom of the plot, up to the corresponding point on
    % the recurrence relation:
    x_coords = [A_0 A_0];
    A_1 = lambda-(K-1)*J^2/A_0;
    y_coords = [A_0 A_1];
    plot(x_coords,y_coords,'k--')
    % This is how you place text labels on the plot:
    %text(mean(x_coords) - 0.035,mean(y_coords),’#1’)
    % Save the current graphics window to a file
    print -depsc cob1.eps

end
%---------------------------------------------------
%ylim([-4, 7]);
%xlim([-5, 5]);
% Put in the third cobweb line. This is a line that goes vertically
% from the previous point up to the recurrence relation.
%x_coords = [P_second P_second];
%P_third = lambda*P_second * (1 - P_second);
%y_coords = [P_second P_third];
%plot(x_coords,y_coords,’b--’)
% Place text label on the plot and save the current graphics window to a file
%text(mean(x_coords) - 0.035,mean(y_coords),’#3’)
%print -depsc cob3.eps

clearvars A_0 A_1 As1 As2 fAs1 fAs2 i J K lambda lim1 lim2 pts steps x_coords y_coords
