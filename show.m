visualizePendulum(t,y(:,1:2))

function visualizePendulum(t, theta)
    figure(3)
    axis equal
    l1=1; l2=1;
    axis([-(l1+l2)*1.2 (l1+l2)*1.2 -(l1+l2)*1.2 (l1+l2)*1.2])
    grid on
    hold on
    title('Pendulum Animation')
    
    % Pendulum model
    x1 = l1 * -sin(theta(:,1));
    y1 = l1 * cos(theta(:,1));
    x2 = x1 + l2 * -sin(theta(:,2));
    y2 = y1 + l2 * cos(theta(:,2));
    pivot = plot(0, 0, 'r.', 'MarkerSize', 20);
    rod1 = line([0 x1(1)], [0 y1(1)], 'LineWidth', 2, 'Color', 'k');
    rod2 = line([x1(1) x2(1)], [y1(1) y2(1)], 'LineWidth', 2, 'Color', 'k');
    bob1 = plot(x1(1), y1(1), 'b.', 'MarkerSize', 20);
    bob2 = plot(x2(1), y2(1), 'g.', 'MarkerSize', 20);

    filename = 'pendulum_animation.gif';
    
    for k = 1:length(t)
        % Update pendulum
        set(rod1, 'XData', [0 x1(k)], 'YData', [0 y1(k)]);
        set(rod2, 'XData', [x1(k) x2(k)], 'YData', [y1(k) y2(k)]);
        set(bob1, 'XData', x1(k), 'YData', y1(k));
        set(bob2, 'XData', x2(k), 'YData', y2(k));

        title(sprintf('Pendulum Animation\nTime = %.2f s, $\\theta_1$ = %.2f rad, $\\theta_2$  %.2f rad', ...
                          t(k), theta(k,1), theta(k,2)), 'Interpreter', 'latex')
        
        drawnow;
        
        frame = getframe(gcf);
        [A, map] = rgb2ind(frame.cdata, 256);
        delay_time = t(end)/length(t); % 10 seconds / 93 frames

        if k == 1
            imwrite(A, map, filename, 'GIF', 'LoopCount', Inf, 'DelayTime', delay_time);
        elseif k~=length(t)
            imwrite(A, map, filename, 'GIF', 'WriteMode', 'append', 'DelayTime', delay_time);
        else
            imwrite(A, map, filename, 'GIF', 'WriteMode', 'append', 'DelayTime', 2); % 2 second final pause
        end
    end
end
