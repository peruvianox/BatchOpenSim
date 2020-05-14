function [Left, Right, Total] = GetMuscleColumns(Data)

% organize columns of muscle outputs from opensim into left and right
% structures


%% Loop through column headers to assign muscles
for i = 1:length(Data.colheaders)
    
    %% initial columns
    if strcmp('time', Data.colheaders(i))  % Time
        Total.Time = Data.data(:,i);
    end
    if strcmp('metabolics_TOTAL', Data.colheaders(i))  % Total Metabolics
        Total.Total = Data.Fdata(:,i);
        Total.Total_UF = Data.data(:,i);
    end
    if strcmp('metabolics_BASAL', Data.colheaders(i)) % Basal Metabolics
        Total.Basal = Data.Fdata(:,i);
        Total.Basal_UF = Data.data(:,i);
    end
    
    %% Left
    % Core muscles
    if strcmp('metabolics_ercspn_l', Data.colheaders(i)) % ErectorSpinae
        Left.ErectorSpinae = Data.Fdata(:,i);
        Left.ErectorSpinae_UF = Data.data(:,i);
        Left.ErectorSpinae_Col = i;
    end
    if strcmp('metabolics_intobl_l', Data.colheaders(i)) % InternalObliques
        Left.InternalObliques = Data.Fdata(:,i);
        Left.InternalObliques_UF = Data.data(:,i);
        Left.InternalObliques_Col = i;
    end
    if strcmp('metabolics_extobl_l', Data.colheaders(i)) % ExternalObliques
        Left.ExternalObliques = Data.Fdata(:,i);
        Left.ExternalObliques_UF = Data.data(:,i);
        Left.ExternalObliques_Col = i;
    end
    
    % Glutes
    if strcmp('metabolics_glut_max1_l', Data.colheaders(i)) % GluteMax1
        Left.GluteMax1= Data.Fdata(:,i);
        Left.GluteMax1_UF = Data.data(:,i);
        Left.GluteMax1_Col = i;
    end
    if strcmp('metabolics_glut_max2_l', Data.colheaders(i)) % GluteMax2
        Left.GluteMax2= Data.Fdata(:,i);
        Left.GluteMax2_UF = Data.data(:,i);
        Left.GluteMax2_Col = i;
    end
    if strcmp('metabolics_glut_max3_l', Data.colheaders(i)) % GluteMax3
        Left.GluteMax3= Data.Fdata(:,i);
        Left.GluteMax3_UF = Data.data(:,i);
        Left.GluteMax3_Col = i;
    end
    if strcmp('metabolics_glut_med1_l', Data.colheaders(i)) % GluteMed1
        Left.GluteMed1= Data.Fdata(:,i);
        Left.GluteMed1_UF = Data.data(:,i);
        Left.GluteMed1_Col = i;
    end
    if strcmp('metabolics_glut_med2_l', Data.colheaders(i)) % GluteMed2
        Left.GluteMed2= Data.Fdata(:,i);
        Left.GluteMed2_UF = Data.data(:,i);
        Left.GluteMed2_Col = i;
    end
    if strcmp('metabolics_glut_med3_l', Data.colheaders(i)) % GluteMed3
        Left.GluteMed3= Data.Fdata(:,i);
        Left.GluteMed3_UF = Data.data(:,i);
        Left.GluteMed3_Col = i;
    end
    if strcmp('metabolics_glut_min1_l', Data.colheaders(i)) % GluteMin1
        Left.GluteMin1= Data.Fdata(:,i);
        Left.GluteMin1_UF = Data.data(:,i);
        Left.GluteMin1_Col = i;
    end
    if strcmp('metabolics_glut_min2_l', Data.colheaders(i)) % GluteMin2
        Left.GluteMin2= Data.Fdata(:,i);
        Left.GluteMin2_UF = Data.data(:,i);
        Left.GluteMin2_Col = i;
    end
    if strcmp('metabolics_glut_min3_l', Data.colheaders(i)) % GluteMin3
        Left.GluteMin3= Data.Fdata(:,i);
        Left.GluteMin3_UF = Data.data(:,i);
        Left.GluteMin3_Col = i;
    end
    
    % Other Hip Muscles
    if strcmp('metabolics_iliacus_l', Data.colheaders(i)) % Iliacus
        Left.Iliacus = Data.Fdata(:,i);
        Left.Iliacus_UF = Data.data(:,i);
        Left.Iliacus_Col = i;
    end
    if strcmp('metabolics_psoas_l', Data.colheaders(i)) % Psoas
        Left.Psoas = Data.Fdata(:,i);
        Left.Psoas_UF = Data.data(:,i);
        Left.Psoas_Col = i;
    end
    if strcmp('metabolics_gem_l', Data.colheaders(i)) % Superior Gemellus
        Left.Gemellus = Data.Fdata(:,i);
        Left.Gemellus_UF = Data.data(:,i);
        Left.Gemellus_Col = i;
    end
    if strcmp('metabolics_peri_l', Data.colheaders(i)) % Piriformis
        Left.Piriformis = Data.Fdata(:,i);
        Left.Piriformis_UF = Data.data(:,i);
        Left.Piriformis_Col = i;
    end
    if strcmp('metabolics_quad_fem_l', Data.colheaders(i)) % Quadratus Femoris
        Left.QuadFem = Data.Fdata(:,i);
        Left.QuadFem_UF = Data.data(:,i);
        Left.QuadFem_Col = i;
    end
    if strcmp('metabolics_tfl_l', Data.colheaders(i)) % Tensor Fascia Latae
        Left.TensorFasciaLatae = Data.Fdata(:,i);
        Left.TensorFasciaLatae_UF = Data.data(:,i);
        Left.TensorFasciaLatae_Col = i;
    end
    
    % Thigh Muscles
    if strcmp('metabolics_semimem_l', Data.colheaders(i)) % Semimembranosis
        Left.Semimem = Data.Fdata(:,i);
        Left.Semimem_UF = Data.data(:,i);
        Left.Semimem_Col = i;
    end
    if strcmp('metabolics_semiten_l', Data.colheaders(i)) % Semitendonosis
        Left.Semiten = Data.Fdata(:,i);
        Left.Semiten_UF = Data.data(:,i);
        Left.Semiten_Col = i;
    end
    if strcmp('metabolics_bifemlh_l', Data.colheaders(i)) % Biceps Femoris Long Head
        Left.BiFemLH = Data.Fdata(:,i);
        Left.BiFemLH_UF = Data.data(:,i);
        Left.BiFemLH_Col = i;
    end
    if strcmp('metabolics_bifemsh_l', Data.colheaders(i)) % Biceps Femoris Short Head
        Left.BiFemSH = Data.Fdata(:,i);
        Left.BiFemSH_UF = Data.data(:,i);
        Left.BiFemSH_Col = i;
    end
    if strcmp('metabolics_sar_l', Data.colheaders(i)) % Sartorius
        Left.Sartorius = Data.Fdata(:,i);
        Left.Sartorius_UF = Data.data(:,i);
        Left.Sartorius_Col = i;
    end
    if strcmp('metabolics_add_long_l', Data.colheaders(i)) % Adductor Longis
        Left.AdductLong = Data.Fdata(:,i);
        Left.AdductLong_UF = Data.data(:,i);
        Left.AdductLong_Col = i;
    end
    if strcmp('metabolics_add_brev_l', Data.colheaders(i)) % Adductor Brevis
        Left.AdductBrev = Data.Fdata(:,i);
        Left.AdductBrev_UF = Data.data(:,i);
        Left.AdductBrev_Col = i;
    end
    if strcmp('metabolics_add_mag1_l', Data.colheaders(i)) % Adductor Magnus 1
        Left.AdductMag1 = Data.Fdata(:,i);
        Left.AdductMag1_UF = Data.data(:,i);
        Left.AdductMag1_Col = i;
    end
    if strcmp('metabolics_add_mag2_l', Data.colheaders(i)) % Adductor Magnus 2
        Left.AdductMag2 = Data.Fdata(:,i);
        Left.AdductMag2_UF = Data.data(:,i);
        Left.AdductMag2_Col = i;
    end
    if strcmp('metabolics_add_mag3_l', Data.colheaders(i)) % Adductor Magnus 3
        Left.AdductMag3 = Data.Fdata(:,i);
        Left.AdductMag3_UF = Data.data(:,i);
        Left.AdductMag3_Col = i;
    end
    if strcmp('metabolics_pect_l', Data.colheaders(i)) % Pectineus
        Left.Pectineus = Data.Fdata(:,i);
        Left.Pectineus_UF = Data.data(:,i);
        Left.Pectineus_Col = i;
    end
    if strcmp('metabolics_grac_l', Data.colheaders(i)) % Gracilis
        Left.Gracilis = Data.Fdata(:,i);
        Left.Gracilis_UF = Data.data(:,i);
        Left.Gracilis_Col = i;
    end
    if strcmp('metabolics_rect_fem_l', Data.colheaders(i)) % Rectus Femoris
        Left.RectFem = Data.Fdata(:,i);
        Left.RectFem_UF = Data.data(:,i);
        Left.RectFem_Col = i;
    end
    if strcmp('metabolics_vas_med_l', Data.colheaders(i)) % Vastus Med
        Left.VastusMed = Data.Fdata(:,i);
        Left.VastusMed_UF = Data.data(:,i);
        Left.VastusMed_Col = i;
    end
    if strcmp('metabolics_vas_int_l', Data.colheaders(i)) % Vastus Int
        Left.VastusInt = Data.Fdata(:,i);
        Left.VastusInt_UF = Data.data(:,i);
        Left.VastusInt_Col = i;
    end
    if strcmp('metabolics_vas_lat_l', Data.colheaders(i)) % Vastus Lat
        Left.VastusLat= Data.Fdata(:,i);
        Left.VastusLat_UF = Data.data(:,i);
        Left.VastusLat_Col = i;
    end
    
    % Shank muscles
    if strcmp('metabolics_soleus_l', Data.colheaders(i)) % Soleus
        Left.Soleus = Data.Fdata(:,i);
        Left.Soleus_UF = Data.data(:,i);
        Left.Soles_Col = i;
    end
    if strcmp('metabolics_med_gas_l', Data.colheaders(i)) % Medial Gastroc
        Left.MedGas = Data.Fdata(:,i);
        Left.MedGas_UF = Data.data(:,i);
        Left.MedGas_Col = i;
    end
    if strcmp('metabolics_lat_gas_l', Data.colheaders(i)) % Lateral Gastroc
        Left.LatGas = Data.Fdata(:,i);
        Left.LatGas_UF = Data.data(:,i);
        Left.LatGas_Col = i;
    end
    if strcmp('metabolics_tib_post_l', Data.colheaders(i)) % tibialis posterior
        Left.TibPost = Data.Fdata(:,i);
        Left.TibPost_UF = Data.data(:,i);
        Left.TibPost_Col = i;
    end
    if strcmp('metabolics_tib_ant_l', Data.colheaders(i)) % tibialis anterior
        Left.TibAnt = Data.Fdata(:,i);
        Left.TibAnt_UF = Data.data(:,i);
        Left.TibAnt_Col = i;
    end
    if strcmp('metabolics_flex_dig_l', Data.colheaders(i)) % Flexor digitorum
        Left.FlexDig = Data.Fdata(:,i);
        Left.FlexDig_UF = Data.data(:,i);
        Left.FlexDig_Col = i;
    end
     if strcmp('metabolics_flex_hal_l', Data.colheaders(i)) % Flexor Hallucius
        Left.FlexHal = Data.Fdata(:,i);
        Left.FlexHal_UF = Data.data(:,i);
        Left.FlexHal_Col = i;
     end
     if strcmp('metabolics_per_brev_l', Data.colheaders(i)) % Peroneus Brevis
         Left.PerBrev = Data.Fdata(:,i);
        Left.PerBrev_UF = Data.data(:,i);
         Left.PerBrev_Col = i;
     end
      if strcmp('metabolics_per_long_l', Data.colheaders(i)) % Peroneus Longus
         Left.PerLong = Data.Fdata(:,i);
        Left.PerLong_UF = Data.data(:,i);
         Left.PerLong_Col = i;
      end
       if strcmp('metabolics_per_tert_l', Data.colheaders(i)) % Peroneus Tert
        Left.PerTert = Data.Fdata(:,i);
        Left.PerTert_UF = Data.data(:,i);
         Left.PerTert_Col = i;
       end
     if strcmp('metabolics_ext_dig_l', Data.colheaders(i)) % Extensor Digitorum
         Left.ExtDig = Data.Fdata(:,i);
         Left.ExtDig_UF = Data.data(:,i);
         Left.ExtDig_Col = i;
     end
     if strcmp('metabolics_ext_hal_l', Data.colheaders(i)) % Extensor Hallucius
         Left.ExtHal = Data.Fdata(:,i);
         Left.ExtHal_UF = Data.data(:,i);
         Left.ExtHal_Col = i;
     end
    
    %% Right
    
     % Core muscles
    if strcmp('metabolics_ercspn_r', Data.colheaders(i)) % ErectorSpinae
        Right.ErectorSpinae = Data.Fdata(:,i);
        Right.ErectorSpinae_UF = Data.data(:,i);
        Right.ErectorSpinae_Col = i;
    end
    if strcmp('metabolics_intobl_r', Data.colheaders(i)) % InternalObliques
        Right.InternalObliques = Data.Fdata(:,i);
        Right.InternalObliques_UF = Data.data(:,i);
        Right.InternalObliques_Col = i;
    end
    if strcmp('metabolics_extobl_r', Data.colheaders(i)) % ExternalObliques
        Right.ExternalObliques = Data.Fdata(:,i);
        Right.ExternalObliques_UF = Data.data(:,i);
        Right.ExternalObliques_Col = i;
    end
    
    % Glutes
    if strcmp('metabolics_glut_max1_r', Data.colheaders(i)) % GluteMax1
        Right.GluteMax1= Data.Fdata(:,i);
        Right.GluteMax1_UF = Data.data(:,i);
        Right.GluteMax1_Col = i;
    end
    if strcmp('metabolics_glut_max2_r', Data.colheaders(i)) % GluteMax2
        Right.GluteMax2= Data.Fdata(:,i);
        Right.GluteMax2_UF = Data.data(:,i);
        Right.GluteMax2_Col = i;
    end
    if strcmp('metabolics_glut_max3_r', Data.colheaders(i)) % GluteMax3
        Right.GluteMax3= Data.Fdata(:,i);
        Right.GluteMax3_UF = Data.data(:,i);
        Right.GluteMax3_Col = i;
    end
    if strcmp('metabolics_glut_med1_r', Data.colheaders(i)) % GluteMed1
        Right.GluteMed1= Data.Fdata(:,i);
        Right.GluteMed1_UF = Data.data(:,i);
        Right.GluteMed1_Col = i;
    end
    if strcmp('metabolics_glut_med2_r', Data.colheaders(i)) % GluteMed2
        Right.GluteMed2= Data.Fdata(:,i);
        Right.GluteMed2_UF = Data.data(:,i);
        Right.GluteMed2_Col = i;
    end
    if strcmp('metabolics_glut_med3_r', Data.colheaders(i)) % GluteMed3
        Right.GluteMed3= Data.Fdata(:,i);
        Right.GluteMed3_UF = Data.data(:,i);
        Right.GluteMed3_Col = i;
    end
    if strcmp('metabolics_glut_min1_r', Data.colheaders(i)) % GluteMin1
        Right.GluteMin1= Data.Fdata(:,i);
        Right.GluteMin1_UF = Data.data(:,i);
        Right.GluteMin1_Col = i;
    end
    if strcmp('metabolics_glut_min2_r', Data.colheaders(i)) % GluteMin2
        Right.GluteMin2= Data.Fdata(:,i);
        Right.GluteMin2_UF = Data.data(:,i);
        Right.GluteMin2_Col = i;
    end
    if strcmp('metabolics_glut_min3_r', Data.colheaders(i)) % GluteMin3
        Right.GluteMin3= Data.Fdata(:,i);
        Right.GluteMin3_UF = Data.data(:,i);
        Right.GluteMin3_Col = i;
    end
    
    % Other Hip Muscles
    if strcmp('metabolics_iliacus_r', Data.colheaders(i)) % Iliacus
        Right.Iliacus = Data.Fdata(:,i);
        Right.Iliacus_UF = Data.data(:,i);
        Right.Iliacus_Col = i;
    end
    if strcmp('metabolics_psoas_r', Data.colheaders(i)) % Psoas
        Right.Psoas = Data.Fdata(:,i);
        Right.Psoas_UF = Data.data(:,i);
        Right.Psoas_Col = i;
    end
    if strcmp('metabolics_gem_r', Data.colheaders(i)) % Superior Gemellus
        Right.Gemellus = Data.Fdata(:,i);
        Right.Gemellus_UF = Data.data(:,i);
        Right.Gemellus_Col = i;
    end
    if strcmp('metabolics_peri_r', Data.colheaders(i)) % Piriformis
        Right.Piriformis = Data.Fdata(:,i);
        Right.Piriformis_UF = Data.data(:,i);
        Right.Piriformis_Col = i;
    end
    if strcmp('metabolics_quad_fem_r', Data.colheaders(i)) % Quadratus Femoris
        Right.QuadFem = Data.Fdata(:,i);
        Right.QuadFem_UF = Data.data(:,i);
        Right.QuadFem_Col = i;
    end
    if strcmp('metabolics_tfl_r', Data.colheaders(i)) % Tensor Fascia Latae
        Right.TensorFasciaLatae = Data.Fdata(:,i);
        Right.TensorFasciaLatae_UF = Data.data(:,i);
        Right.TensorFasciaLatae_Col = i;
    end
    
    % Thigh Muscles
    if strcmp('metabolics_semimem_r', Data.colheaders(i)) % Semimembranosis
        Right.Semimem = Data.Fdata(:,i);
        Right.Semimem_UF = Data.data(:,i);
        Right.Semimem_Col = i;
    end
    if strcmp('metabolics_semiten_r', Data.colheaders(i)) % Semitendonosis
        Right.Semiten = Data.Fdata(:,i);
        Right.Semiten_UF = Data.data(:,i);
        Right.Semiten_Col = i;
    end
    if strcmp('metabolics_bifemlh_r', Data.colheaders(i)) % Biceps Femoris Long Head
        Right.BiFemLH = Data.Fdata(:,i);
        Right.BiFemLH_UF = Data.data(:,i);
        Right.BiFemLH_Col = i;
    end
    if strcmp('metabolics_bifemsh_r', Data.colheaders(i)) % Biceps Femoris Short Head
        Right.BiFemSH = Data.Fdata(:,i);
        Right.BiFemSH_UF = Data.data(:,i);
        Right.BiFemSH_Col = i;
    end
    if strcmp('metabolics_sar_r', Data.colheaders(i)) % Sartorius
        Right.Sartorius = Data.Fdata(:,i);
        Right.Sartorius_UF = Data.data(:,i);
        Right.Sartorius_Col = i;
    end
    if strcmp('metabolics_add_long_r', Data.colheaders(i)) % Adductor Longis
        Right.AdductLong = Data.Fdata(:,i);
        Right.AdductLong_UF = Data.data(:,i);
        Right.AdductLong_Col = i;
    end
    if strcmp('metabolics_add_brev_r', Data.colheaders(i)) % Adductor Brevis
        Right.AdductBrev = Data.Fdata(:,i);
        Right.AdductBrev_UF = Data.data(:,i);
        Right.AdductBrev_Col = i;
    end
    if strcmp('metabolics_add_mag1_r', Data.colheaders(i)) % Adductor Magnus 1
        Right.AdductMag1 = Data.Fdata(:,i);
        Right.AdductMag1_UF = Data.data(:,i);
        Right.AdductMag1_Col = i;
    end
    if strcmp('metabolics_add_mag2_r', Data.colheaders(i)) % Adductor Magnus 2
        Right.AdductMag2 = Data.Fdata(:,i);
        Right.AdductMag2_UF = Data.data(:,i);
        Right.AdductMag2_Col = i;
    end
    if strcmp('metabolics_add_mag3_r', Data.colheaders(i)) % Adductor Magnus 3
        Right.AdductMag3 = Data.Fdata(:,i);
        Right.AdductMag3_UF = Data.data(:,i);
        Right.AdductMag3_Col = i;
    end
    if strcmp('metabolics_pect_r', Data.colheaders(i)) % Pectineus
        Right.Pectineus = Data.Fdata(:,i);
        Right.Pectineus_UF = Data.data(:,i);
        Right.Pectineus_Col = i;
    end
    if strcmp('metabolics_grac_r', Data.colheaders(i)) % Gracilis
        Right.Gracilis = Data.Fdata(:,i);
        Right.Gracilis_UF = Data.data(:,i);
        Right.Gracilis_Col = i;
    end
    if strcmp('metabolics_rect_fem_r', Data.colheaders(i)) % Rectus Femoris
        Right.RectFem = Data.Fdata(:,i);
        Right.RectFem_UF = Data.data(:,i);
        Right.RectFem_Col = i;
    end
    if strcmp('metabolics_vas_med_r', Data.colheaders(i)) % Vastus Med
        Right.VastusMed = Data.Fdata(:,i);
        Right.VastusMed_UF = Data.data(:,i);
        Right.VastusMed_Col = i;
    end
    if strcmp('metabolics_vas_int_r', Data.colheaders(i)) % Vastus Int
        Right.VastusInt = Data.Fdata(:,i);
        Right.VastusInt_UF = Data.data(:,i);
        Right.VastusInt_Col = i;
    end
    if strcmp('metabolics_vas_lat_r', Data.colheaders(i)) % Vastus Lat
        Right.VastusLat= Data.Fdata(:,i);
        Right.VastusLat_UF = Data.data(:,i);
        Right.VastusLat_Col = i;
    end
    
    % Shank muscles
    if strcmp('metabolics_soleus_r', Data.colheaders(i)) % Soleus
        Right.Soleus = Data.Fdata(:,i);
        Right.Soleus_UF = Data.data(:,i);
        Right.Soles_Col = i;
    end
    if strcmp('metabolics_med_gas_r', Data.colheaders(i)) % Medial Gastroc
        Right.MedGas = Data.Fdata(:,i);
        Right.MedGas_UF = Data.data(:,i);
        Right.MedGas_Col = i;
    end
    if strcmp('metabolics_lat_gas_r', Data.colheaders(i)) % Lateral Gastroc
        Right.LatGas = Data.Fdata(:,i);
        Right.LatGas_UF = Data.data(:,i);
        Right.LatGas_Col = i;
    end
    if strcmp('metabolics_tib_post_r', Data.colheaders(i)) % tibialis posterior
        Right.TibPost = Data.Fdata(:,i);
        Right.TibPost_UF = Data.data(:,i);
        Right.TibPost_Col = i;
    end
    if strcmp('metabolics_tib_ant_r', Data.colheaders(i)) % tibialis anterior
        Right.TibAnt = Data.Fdata(:,i);
        Right.TibAnt_UF = Data.data(:,i);
        Right.TibAnt_Col = i;
    end
    if strcmp('metabolics_flex_dig_r', Data.colheaders(i)) % Flexor digitorum
        Right.FlexDig = Data.Fdata(:,i);
        Right.FlexDig_UF = Data.data(:,i);
        Right.FlexDig_Col = i;
    end
     if strcmp('metabolics_flex_hal_r', Data.colheaders(i)) % Flexor Hallucius
        Right.FlexHal = Data.Fdata(:,i);
        Right.FlexHal_UF = Data.data(:,i);
        Right.FlexHal_Col = i;
     end
     if strcmp('metabolics_per_brev_r', Data.colheaders(i)) % Peroneus Brevis
         Right.PerBrev = Data.Fdata(:,i);
        Right.PerBrev_UF = Data.data(:,i);
         Right.PerBrev_Col = i;
     end
      if strcmp('metabolics_per_long_r', Data.colheaders(i)) % Peroneus Longus
         Right.PerLong = Data.Fdata(:,i);
        Right.PerLong_UF = Data.data(:,i);
         Right.PerLong_Col = i;
      end
       if strcmp('metabolics_per_tert_r', Data.colheaders(i)) % Peroneus Tert
        Right.PerTert = Data.Fdata(:,i);
        Right.PerTert_UF = Data.data(:,i);
         Right.PerTert_Col = i;
       end
     if strcmp('metabolics_ext_dig_r', Data.colheaders(i)) % Extensor Digitorum
         Right.ExtDig = Data.Fdata(:,i);
         Right.ExtDig_UF = Data.data(:,i);
         Right.ExtDig_Col = i;
     end
     if strcmp('metabolics_ext_hal_r', Data.colheaders(i)) % Extensor Hallucius
         Right.ExtHal = Data.Fdata(:,i);
         Right.ExtHal_UF = Data.data(:,i);
         Right.ExtHal_Col = i;
     end
    
    
    
%     % Soleus
%     if strcmp('metabolics_soleus_r', Data.colheaders(i))
%         Right.Soleus = Data.Fdata(:,i);
%         Right.Soleus_UF = Data.data(:,i);
%     end
%     % Medial Gastroc
%     if strcmp('metabolics_med_gas_r', Data.colheaders(i))
%         Right.MedGas = Data.Fdata(:,i);
%         Right.MedGas_UF = Data.data(:,i);
%     end
%     
%     % Lateral Gastroc
%     if strcmp('metabolics_lat_gas_r', Data.colheaders(i))
%         Right.LatGas = Data.Fdata(:,i);
%         Right.LatGas_UF = Data.data(:,i);
%     end
    
    
    
end


