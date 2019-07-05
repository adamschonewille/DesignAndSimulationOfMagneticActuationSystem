function [Br_min Br_nom] = magnetLookup( grade )
% AUTHOR:   Adam Schonewille
% DATE:     October 17 2018
% ABOUT:    Finds the nominal and maximum Remanence values of a NdFeB
%       magnet from a lookup table
%
%% INPUTS
% grade  - [string] The of magnet that is being designed for
%
%% RETURNED OUTPUTS
% Br_min - [Teslas] The minimum Remanence value for the magnet 
% Br_nom - [Teslas] The nominal Remanence value for the magnet

%% Function Code
% Table created from Eclipse Magnetics Tables on:
% Sintered Neodymium Iron Boron (NdFeB) Magnets [mT]
gradeTable = [  "N35"   1170 1220;
                "N38"   1220 1260;
                "N40"   1260 1300;
                "N42"   1300 1330;
                "N45"   1330 1370;
                "N48"   1370 1410;
                "N50"   1410 1440;
                "N52"   1440 1470;
                "N33M"  1140 1170;
                "N35M"  1170 1220;
                "N38M"  1220 1260;
                "N40M"  1260 1300;
                "N42M"  1300 1330;
                "N45M"  1330 1370;
                "N48M"  1370 1410;
                "N50M"  1410 1440;
                "N30H"  1080 1140;
                "N33H"  1140 1170;
                "N35H"  1170 1220;
                "N38H"  1220 1260;
                "N40H"  1260 1300;
                "N42H"  1300 1330;
                "N44H"  1330 1360;
                "N46H"  1360 1380;
                "N48H"  1380 1410;
                "N30SH" 1080 1140;
                "N33SH" 1140 1170;
                "N35SH" 1170 1220;
                "N38SH" 1220 1260;
                "N40SH" 1260 1300;
                "N42SH" 1300 1330;
                "N44SH" 1330 1360;
                "N46SH" 1360 1380;
                "N28UH" 1040 1080;
                "N30UH" 1080 1140;
                "N33UH" 1140 1170;
                "N35UH" 1170 1220;
                "N38UH" 1220 1260;
                "N40UH" 1260 1300;
                "N28EH" 1040 1080;
                "N30EH" 1080 1140;
                "N33EH" 1140 1170;
                "N35EH" 1170 1220;
                "N38EH" 1220 1260;
                "N25AH" 970  1020];
A = find(gradeTable==grade);
    if ( isempty(A) ) 
        error("No valid magnet grade requested");
    elseif ( A > length(gradeTable) )
        error("No valid magnet grade requested");    
    else
        Br_min = gradeTable(A,2);
        Br_min = str2num(Br_min)/1000;
        Br_nom = gradeTable(A,3);
        Br_nom = str2num(Br_nom)/1000;
    end
    
end