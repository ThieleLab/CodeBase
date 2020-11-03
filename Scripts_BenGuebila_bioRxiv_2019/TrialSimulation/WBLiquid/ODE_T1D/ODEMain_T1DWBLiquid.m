function [tout, yout] = ODEMain_T1DWBLiquid

    tic

    clear global;

    global eventsFunctionHandle;
    eventsFunctionHandle = @events;

    SetupTableParameters;

    y0 = ODEInitialValues;
    tstart = 0;
    tout = tstart;

    %call ODE RHS function in order to init all parameters
    ODERHSFunction(tstart, y0);
    %perform initial switches update
    [y0 switchUpdate] = PerformSwitches(tstart, y0);

    yout = y0.';

%     outtimes = [ 0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360 375 390 405 420 435 450 465 480 495 510 525 540 555 570 585 600 615 630 645 660 675 690 705 720 735 750 765 780 795 810 825 840 855 870 885 900 915 930 945 960 975 990 1005 1020 1035 1050 1065 1080 1095 1110 1125 1140 1155 1170 1185 1200 1215 1230 1245 1260 1275 1290 1305 1320 1335 1350 1365 1380 1395 1410 1425 1440 ];
    outtimes = [0:100:900 1000+[0:5:600]];
        
    global switchtimes;
    switchtimes = [ 0 15 180 1000 ];

    NextRestartTime = tstart;
    while NextRestartTime < outtimes(end)

        switchtimes = switchtimes(switchtimes > NextRestartTime);
        if isempty(switchtimes)
            NextRestartTime = outtimes(end);
        else
            NextRestartTime = switchtimes(1);
        end

        while (1)
            outtimes_step_i = [tstart outtimes(outtimes > tstart & outtimes <= NextRestartTime)];
            if (outtimes_step_i(end) < NextRestartTime)
                outtimes_step_i = [outtimes_step_i NextRestartTime]; %#ok<AGROW>
            end

            [t, y] = ode15s(@ODERHSFunction, outtimes_step_i, y0, ODEoptions);

            nt = length(t);
            tout = [tout; t(2:nt)]; %#ok<AGROW>
            yout = [yout; y(2:nt,:)]; %#ok<AGROW>

            y0 = y(nt,:);
            tstart = t(nt);

            %update start vector by switches (if applies)
            [y0 switchUpdate] = PerformSwitches(tstart, y0);

            if (tstart==outtimes_step_i(end))
                break
            end
        end
    end

    [tout idx_out] = intersect(tout, outtimes);
    yout = yout(idx_out, :);

    toc

function [value,isterminal,direction] = events(Time,y)

    value = 1;      %no event per default
    isterminal = 1; % stop the integration
    direction = 0;  % negative direction

    global switchtimes;
    if (min(abs(switchtimes - Time)) == 0.0)
        return
    end

    [yOut switchUpdate] = PerformSwitches(Time,y);

    if ~switchUpdate
        %no switch fired
        return
    end

    value = 0;      % force system restart
