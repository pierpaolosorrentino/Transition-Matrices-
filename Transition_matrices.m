%% load avalanches
% insert path to avalanches. Avalanches are stored in a cell array
% containing n cells, where n is the number of subjects. Each cell
% contains, in turn, one more cell array with m cells, with m = the number
% of individual avalanches. Each cell contains a matrix with raws as
% regions of interest, and columns as timesteps of the specific avalanche.

avalanches=load('');

%% load tractography
% dti's are loaded as a 3D matrix, where raws and columns are regions of
% interest, while subjects are along the third dimension.

DTIs=load('');
dti_media=mean(DTIs,3);

%% parameters

num_reg=90; % number of regions of interest 
group_generator = @(x)x(randperm(numel(x))); % anonimous function to generate random positions
num_perm=10; % number of random surrogates that will be generated
min_size_aval=30; % selects the minimum avalanche size to be taken into account

%%
mask = tril(true(size(dti_media)),-1);  % to masks matrices 
dti_media_vect = dti_media(mask); % to mask matrices

% computation of transition matrices
sum_per_group=zeros(num_reg,num_reg);
    for kk1=1:size(avalanches,2)
        sum_per_pat=zeros(num_reg,num_reg);
        aval_longer_than_min=1;
        for kk2=1:size(avalanches{kk1},2)
            if size(avalanches{kk1}{kk2},2)>min_size_aval % select minimal size of avalanches
            curr_aval=avalanches{kk1}{kk2};
            temp=func_transition_matrix(curr_aval);
            sum_per_pat=sum_per_pat+temp;
            aval_longer_than_min=aval_longer_than_min+1;
            end    
        end
        prob_per_pat(:,:,kk1)=sum_per_pat./aval_longer_than_min; % normalizes by the number of avalanches
    end

P1=mean(prob_per_pat,3);
P1=P1+P1'./2; % symmetrizes the average transition matrix
BB = P1(mask);  
clearvars temp_pos

[r,p]=corr(dti_media_vect,BB,'Type','Spearman');

% generates surrogate correlations by randomizing the transition matrix
    for kk1=1:num_perm
         temp_pos=group_generator(1:size(BB,1));   
         temp=BB(temp_pos);  
         [r_rand(kk1),p_rand(kk1)]= corr(dti_media_vect,temp,'Type','Spearman');
    end
    
% corrected p's and r's    
p_corr=length(find(p>p_rand))/num_perm;
% r_corr=length(find(abs(r)<abs(r_rand)))/num_perm;
    
clearvars -except prob_per_pat avalanches group_generator dti_media num_perm mask r_corr p_corr dti_media_vect min_size_aval r p p_rand r_rand P1 num_reg
     
% computation of randomized matrices
sum_per_goup=zeros(num_reg,num_reg);

for kk1=1:size(avalanches,2)
    sum_per_pat=zeros(num_reg,num_reg);
    c=cellfun('size',avalanches{kk1},2);
    long_aval=sum(c>min_size_aval);
    if long_aval>0
    big_avals=1;
    for kk2=1:size(avalanches{kk1},2)
        if size(avalanches{kk1}{kk2},2)>min_size_aval % selection of minimal avalanche
            curr_aval=avalanches{kk1}{kk2};
            temp=zeros(size(curr_aval,1),size(curr_aval,1),num_perm);
            for kk3=1:num_perm
                rand_temp=group_generator(1:size(curr_aval,2));
                curr_aval_rand=curr_aval(:,rand_temp); % randomizes temporal order or avalanches
                temp(:,:,kk3)=func_transition_matrix(curr_aval_rand);
            end
            tm_per_aval_t(:,:,big_avals)=reshape(permute(temp,[1 3 2]),[],size(temp,2),1);
            big_avals=big_avals+1;
        end
    end  
    tm_per_pat_t(:,:,kk1)=mean(tm_per_aval_t,3);% block matrix=raws as num regs x numb of surrogates, cols= numb or regions, 3rd dimension is numb of subjects
    clearvars  hh_struct tm_per_aval_t   
    end
    disp(kk1) % display number of subject
end

tm_global=mean(tm_per_pat_t,3);
hh_struc=mat2cell(tm_global,repelem(90,num_perm),90); 
a=cat(3,hh_struc{:}); % stores surrogate matrices as region x region x number of surrogate (permutation). i.e. each random surroate is the average across subjects

for kk1=1:size(a,3)
   probe = a(:,:,kk1); 
   probe=probe+probe'./2;
   BB=probe(mask);
   [r_rand_corr(kk1),p_rand_corr(kk1)]=corr(dti_media_vect,BB,'Type','Spearman');
end  

final_p=length(find(p_corr>p_rand_corr))/num_perm;
% final_r=length(find(abs(r)<abs(r_rand_corr)))/num_perm;
