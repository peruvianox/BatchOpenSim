function[Energy] = Integrate4EnergyCost(Left, Right, LeftTot, RightTot, LogTimes_L, LogTimes_R, WalkSpeed)


SampFreq = 1000; 
LStride = sum(LogTimes_L.Stride); 
RStride = sum(LogTimes_R.Stride); 
% LStride = sum(LogTimes_L.Stride) / SampFreq; 
% RStride = sum(LogTimes_R.Stride) / SampFreq; 
LStrideLength = LStride * WalkSpeed;
RStrideLength = RStride * WalkSpeed;

% Full Stride Integrals and time calculations
Energy.Stride.Left_Time = length(Left(LogTimes_L.Stride)) ./ SampFreq; 
% Energy.Stride.Left_Muscles_Int = trapz(Left(LogTimes_L.Stride, :)) / SampFreq / LStride;
Energy.Stride.Left_Muscles_Int = trapz(Left(LogTimes_L.Stride, :)) / LStride;
% Energy.Stride.Left_Muscles_Int = trapz(Left(LogTimes_L.Stride, :)) / LStrideLength;
Energy.Stride.Right_Time = length(Right(LogTimes_R.Stride)) ./ SampFreq;
% Energy.Stride.Right_Muscles_Int = trapz(Right(LogTimes_R.Stride, :)) / RStride / SampFreq;
Energy.Stride.Right_Muscles_Int = trapz(Right(LogTimes_R.Stride, :)) / RStride;
% Energy.Stride.Right_Muscles_Int = trapz(Right(LogTimes_R.Stride, :)) / RStrideLength;
Energy.Stride.Avg_Muscles_Int = mean([Energy.Stride.Left_Muscles_Int; Energy.Stride.Right_Muscles_Int]);

% Energy.Stride.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stride)) / SampFreq / LStride;
% Energy.Stride.Right_Total_Int = trapz(RightTot(LogTimes_R.Stride)) / RStride / SampFreq;
Energy.Stride.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stride)) / LStride;
Energy.Stride.Right_Total_Int = trapz(RightTot(LogTimes_R.Stride)) / RStride;
% Energy.Stride.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stride)) / LStrideLength;
% Energy.Stride.Right_Total_Int = trapz(RightTot(LogTimes_R.Stride)) / RStrideLength;
Energy.Stride.Avg_Total_Int = mean([Energy.Stride.Left_Total_Int; Energy.Stride.Right_Total_Int]); 

% integreate for specific stride durations
% stance
Energy.Stance.Left_Time = length(Left(LogTimes_L.Stance)) ./ SampFreq; 
% Energy.Stance.Left_Muscles_Int = trapz(Left(LogTimes_L.Stance, :)) / SampFreq / LStride;
Energy.Stance.Left_Muscles_Int = trapz(Left(LogTimes_L.Stance, :)) / LStride;
% Energy.Stance.Left_Muscles_Int = trapz(Left(LogTimes_L.Stance, :)) / LStrideLength;
Energy.Stance.Right_Time = length(Right(LogTimes_R.Stance)) ./ SampFreq;
% Energy.Stance.Right_Muscles_Int = trapz(Right(LogTimes_R.Stance, :)) / RStride / SampFreq;
Energy.Stance.Right_Muscles_Int = trapz(Right(LogTimes_R.Stance, :)) / RStride;
% Energy.Stance.Right_Muscles_Int = trapz(Right(LogTimes_R.Stance, :)) / RStrideLength;
Energy.Stance.Avg_Muscles_Int = mean([Energy.Stance.Left_Muscles_Int; Energy.Stance.Right_Muscles_Int]);

% Energy.Stance.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stance)) / SampFreq / LStride;
% Energy.Stance.Right_Total_Int = trapz(RightTot(LogTimes_R.Stance)) / RStride / SampFreq;
Energy.Stance.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stance)) / LStride;
Energy.Stance.Right_Total_Int = trapz(RightTot(LogTimes_R.Stance)) / RStride;
% Energy.Stance.Left_Total_Int = trapz(LeftTot(LogTimes_L.Stance)) / LStrideLength;
% Energy.Stance.Right_Total_Int = trapz(RightTot(LogTimes_R.Stance)) / RStrideLength;
Energy.Stance.Avg_Total_Int = mean([Energy.Stance.Left_Total_Int; Energy.Stance.Right_Total_Int]); 

% Initial double support
Energy.DS1.Left_Time = length(Left(LogTimes_L.DS1)) ./ SampFreq; 
% Energy.DS1.Left_Muscles_Int = trapz(Left(LogTimes_L.DS1, :)) / SampFreq / LStride;
Energy.DS1.Left_Muscles_Int = trapz(Left(LogTimes_L.DS1, :)) / LStride;
% Energy.DS1.Left_Muscles_Int = trapz(Left(LogTimes_L.DS1, :)) / LStrideLength;
Energy.DS1.Right_Time = length(Right(LogTimes_R.DS1)) ./ SampFreq;
% Energy.DS1.Right_Muscles_Int = trapz(Right(LogTimes_R.DS1, :)) / RStride / SampFreq;
Energy.DS1.Right_Muscles_Int = trapz(Right(LogTimes_R.DS1, :)) / RStride;
% Energy.DS1.Right_Muscles_Int = trapz(Right(LogTimes_R.DS1, :)) / RStrideLength;
Energy.DS1.Avg_Muscles_Int = mean([Energy.DS1.Left_Muscles_Int; Energy.DS1.Right_Muscles_Int]);

% Energy.DS1.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS1)) / SampFreq / LStride;
% Energy.DS1.Right_Total_Int = trapz(RightTot(LogTimes_R.DS1)) / RStride / SampFreq;
Energy.DS1.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS1)) / LStride;
Energy.DS1.Right_Total_Int = trapz(RightTot(LogTimes_R.DS1)) / RStride;
% Energy.DS1.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS1)) / LStrideLength;
% Energy.DS1.Right_Total_Int = trapz(RightTot(LogTimes_R.DS1)) / RStrideLength;
Energy.DS1.Avg_Total_Int = mean([Energy.DS1.Left_Total_Int; Energy.DS1.Right_Total_Int]); 

% single support
Energy.SingSup.Left_Time = length(Left(LogTimes_L.SingSup)) ./ SampFreq; 
% Energy.SingSup.Left_Muscles_Int = trapz(Left(LogTimes_L.SingSup, :)) / SampFreq / LStride;
Energy.SingSup.Left_Muscles_Int = trapz(Left(LogTimes_L.SingSup, :)) / LStride;
% Energy.SingSup.Left_Muscles_Int = trapz(Left(LogTimes_L.SingSup, :)) / LStrideLength;
Energy.SingSup.Right_Time = length(Right(LogTimes_R.SingSup)) ./ SampFreq;
% Energy.SingSup.Right_Muscles_Int = trapz(Right(LogTimes_R.SingSup, :)) / RStride / SampFreq;
Energy.SingSup.Right_Muscles_Int = trapz(Right(LogTimes_R.SingSup, :)) / RStride;
% Energy.SingSup.Right_Muscles_Int = trapz(Right(LogTimes_R.SingSup, :)) / RStrideLength;
Energy.SingSup.Avg_Muscles_Int = mean([Energy.SingSup.Left_Muscles_Int; Energy.SingSup.Right_Muscles_Int]);

% Energy.SingSup.Left_Total_Int = trapz(LeftTot(LogTimes_L.SingSup)) / SampFreq / LStride;
% Energy.SingSup.Right_Total_Int = trapz(RightTot(LogTimes_R.SingSup)) / RStride / SampFreq;
Energy.SingSup.Left_Total_Int = trapz(LeftTot(LogTimes_L.SingSup)) / LStride;
Energy.SingSup.Right_Total_Int = trapz(RightTot(LogTimes_R.SingSup)) / RStride;
% Energy.SingSup.Left_Total_Int = trapz(LeftTot(LogTimes_L.SingSup)) / LStrideLength;
% Energy.SingSup.Right_Total_Int = trapz(RightTot(LogTimes_R.SingSup)) / RStrideLength;
Energy.SingSup.Avg_Total_Int = mean([Energy.SingSup.Left_Total_Int; Energy.SingSup.Right_Total_Int]); 

% final double support
Energy.DS2.Left_Time = length(Left(LogTimes_L.DS2)) ./ SampFreq; 
% Energy.DS2.Left_Muscles_Int = trapz(Left(LogTimes_L.DS2, :)) / SampFreq / LStride;
Energy.DS2.Left_Muscles_Int = trapz(Left(LogTimes_L.DS2, :)) / LStride;
% Energy.DS2.Left_Muscles_Int = trapz(Left(LogTimes_L.DS2, :)) / LStrideLength;
Energy.DS2.Right_Time = length(Right(LogTimes_R.DS2)) ./ SampFreq;
% Energy.DS2.Right_Muscles_Int = trapz(Right(LogTimes_R.DS2, :)) / RStride / SampFreq;
Energy.DS2.Right_Muscles_Int = trapz(Right(LogTimes_R.DS2, :)) / RStride;
% Energy.DS2.Right_Muscles_Int = trapz(Right(LogTimes_R.DS2, :)) / RStrideLength;
Energy.DS2.Avg_Muscles_Int = mean([Energy.DS2.Left_Muscles_Int; Energy.DS2.Right_Muscles_Int]);

% Energy.DS2.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS2)) / SampFreq / LStride;
% Energy.DS2.Right_Total_Int = trapz(RightTot(LogTimes_R.DS2)) / RStride / SampFreq;
Energy.DS2.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS2)) / LStride;
Energy.DS2.Right_Total_Int = trapz(RightTot(LogTimes_R.DS2)) / RStride;
% Energy.DS2.Left_Total_Int = trapz(LeftTot(LogTimes_L.DS2)) / LStrideLength;
% Energy.DS2.Right_Total_Int = trapz(RightTot(LogTimes_R.DS2)) / RStrideLength;
Energy.DS2.Avg_Total_Int = mean([Energy.DS2.Left_Total_Int; Energy.DS2.Right_Total_Int]); 

% swing phase
Energy.Swing.Left_Time = length(Left(LogTimes_L.Swing)) ./ SampFreq; 
% Energy.Swing.Left_Muscles_Int = trapz(Left(LogTimes_L.Swing, :)) / SampFreq / LStride;
Energy.Swing.Left_Muscles_Int = trapz(Left(LogTimes_L.Swing, :)) / LStride;
% Energy.Swing.Left_Muscles_Int = trapz(Left(LogTimes_L.Swing, :)) / LStrideLength;
Energy.Swing.Right_Time = length(Right(LogTimes_R.Swing)) ./ SampFreq;
% Energy.Swing.Right_Muscles_Int = trapz(Right(LogTimes_R.Swing, :)) / RStride / SampFreq;
Energy.Swing.Right_Muscles_Int = trapz(Right(LogTimes_R.Swing, :)) / RStride;
% Energy.Swing.Right_Muscles_Int = trapz(Right(LogTimes_R.Swing, :)) / RStrideLength;
Energy.Swing.Avg_Muscles_Int = mean([Energy.Swing.Left_Muscles_Int; Energy.Swing.Right_Muscles_Int]);


% Energy.Swing.Right_Total_Int = trapz(RightTot(LogTimes_R.Swing)) / RStride / SampFreq;
% Energy.Swing.Left_Total_Int = trapz(LeftTot(LogTimes_L.Swing)) / SampFreq / LStride;
Energy.Swing.Left_Total_Int = trapz(LeftTot(LogTimes_L.Swing)) / LStride;
Energy.Swing.Right_Total_Int = trapz(RightTot(LogTimes_R.Swing)) / RStride;
% Energy.Swing.Left_Total_Int = trapz(LeftTot(LogTimes_L.Swing)) / LStrideLength;
% Energy.Swing.Right_Total_Int = trapz(RightTot(LogTimes_R.Swing)) / RStrideLength;
Energy.Swing.Avg_Total_Int = mean([Energy.Swing.Left_Total_Int; Energy.Swing.Right_Total_Int]); 

end



