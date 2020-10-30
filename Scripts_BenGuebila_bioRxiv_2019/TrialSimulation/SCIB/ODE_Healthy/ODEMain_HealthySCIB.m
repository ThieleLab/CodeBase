function [tout, yout] = ODEMain_HealthySCIB

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

%     outtimes = [ 0 7.5 15 22.5 30 37.5 45 52.5 60 67.5 75 82.5 90 97.5 105 112.5 120 127.5 135 142.5 150 157.5 165 172.5 180 187.5 195 202.5 210 217.5 225 232.5 240 247.5 255 262.5 270 277.5 285 292.5 300 307.5 315 322.5 330 337.5 345 352.5 360 367.5 375 382.5 390 397.5 405 412.5 420 427.5 435 442.5 450 457.5 465 472.5 480 487.5 495 502.5 510 517.5 525 532.5 540 547.5 555 562.5 570 577.5 585 592.5 600 607.5 615 622.5 630 637.5 645 652.5 660 667.5 675 682.5 690 697.5 705 712.5 720 727.5 735 742.5 750 757.5 765 772.5 780 787.5 795 802.5 810 817.5 825 832.5 840 847.5 855 862.5 870 877.5 885 892.5 900 907.5 915 922.5 930 937.5 945 952.5 960 967.5 975 982.5 990 997.5 1005 1012.5 1020 1027.5 1035 1042.5 1050 1057.5 1065 1072.5 1080 1087.5 1095 1102.5 1110 1117.5 1125 1132.5 1140 1147.5 1155 1162.5 1170 1177.5 1185 1192.5 1200 1207.5 1215 1222.5 1230 1237.5 1245 1252.5 1260 1267.5 1275 1282.5 1290 1297.5 1305 1312.5 1320 1327.5 1335 1342.5 1350 1357.5 1365 1372.5 1380 1387.5 1395 1402.5 1410 1417.5 1425 1432.5 1440 ];
    outtimes = [0:100:900 1000+[0:5:600]];
        
    global switchtimes;
    switchtimes = [ 0 15 180 1000 1015 ];

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
