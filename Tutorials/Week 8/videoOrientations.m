function videoOrientations(AHRS,videoOutputFilename)

telapsed = 0;

h = figure;
nextFrameTime = 0; % s
frameRate = 30; % fps

for i = 1:(size(AHRS.R,3)-1)

    nextdt = AHRS.t(i+1)-AHRS.t(i);
    
    % Fast moving orientation
    plot3([0 AHRS.R(1,1,i)],[0 AHRS.R(2,1,i)],[0 AHRS.R(3,1,i)],'r','LineWidth',3);% Acc commponents
    hold on;
    plot3([0 AHRS.R(1,2,i)],[0 AHRS.R(2,2,i)],[0 AHRS.R(3,2,i)],'g','LineWidth',3);% Orth mag commponents
    plot3([0 AHRS.R(1,3,i)],[0 AHRS.R(2,3,i)],[0 AHRS.R(3,3,i)],'b','LineWidth',3);% third basis vector

    scaleG = 1/norm(AHRS.Acc(i,:));
    scaleM = 1/norm(AHRS.Mag(i,:));
    plot3([0 scaleG*AHRS.Acc(i,1)],[0 scaleG*AHRS.Acc(i,2)],[0 scaleG*AHRS.Acc(i,3)],'y','LineWidth',4);% third basis vector
    plot3([0 scaleM*AHRS.Mag(i,1)],[0 scaleM*AHRS.Mag(i,2)],[0 0],'m','LineWidth',4);% third basis vector
    
    scalePatchX = 0.25; % Length of phone in x direction
    scalePatchY = 0.75; % Length of phone in y direction
    posXvec = scalePatchX*AHRS.R(:,1,i);
    posYvec = scalePatchY*AHRS.R(:,2,i);
    patchVectors = [posXvec+posYvec , posXvec-posYvec, -posXvec-posYvec, -posXvec+posYvec ]';
    patch(patchVectors(:,1),patchVectors(:,2),patchVectors(:,3),[0.5 0.5 0.5],'FaceAlpha',0.7)
    
    hold off;
    
    xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);    
    az = 30;
%     az = mod(0.01*360*telapsed,360); 
    el = 15; 
    view(az,el);
    xlabel('Global x');ylabel('Global y');zlabel('Global z ');
    axis('square');
    grid on;
    legend({'Device x (out of right side)','Device y (up along glass)','Device z (out of face)','Estimated acceleration (normalised)','Estiamted xy of magnetic field (normalised)'})
    title(['Time (s): ',num2str(AHRS.t(i))],'Position',[1 1 1.05])

    drawnow;
    
    if i == 1
        disp('Prepare figure for capture (resize, etc.), then press any key to continue...')
        pause;
    end
    
    % Great the video and/or add another frame to it
    if i == 1;
        v = VideoWriter(videoOutputFilename);
        open(v);
    end
    if AHRS.t(i) >= nextFrameTime % Wait for time when next frame starts
        numFramesSkipped = floor(AHRS.t(i) - nextFrameTime) / (1/frameRate);
        if numFramesSkipped > 0 % No frames skipped
            for k = 1:numFramesSkipped % Do skipped frames
                writeVideo(v,oldframe); % Fill with old frames
            end
        end
        frame = getframe(h); % get most recent frame
        writeVideo(v,frame);
        nextFrameTime = nextFrameTime + (1/frameRate);
        oldframe = frame;
    end

    
    telapsed = telapsed + nextdt;
    
end
close(v);
