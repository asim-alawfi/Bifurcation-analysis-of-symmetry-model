function grid=dde_coll_set_grid(purpose,degree,varargin)
%% define various node grids for polynomial interpolation
% $Id: dde_coll_set_grid.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'type',[],'subset',@(x)x(2:end),'force_smooth',false,...
    'lhs_num',1,'daetest','rank','maxdegree_linear',10};
options=dde_set_options(default,varargin,'pass_on');
if isa(options.type,'function_handle')
    grid=options.type(purpose,degree);
    return
elseif isnumeric(options.type) && ~isempty(options.type)
    grid=options.type;
    return
end
dim=size(options.lhs_num,1);
switch options.daetest
    case 'rank'
        isdae=rank(options.lhs_num)<dim;
    case 'false'
        isdae=false;
    otherwise
        isdae=any(options.lhs_num(:)~=reshape(eye(dim),[],1));
end
if isempty(options.type) % ||strcmp(options.type,'force_smooth')
    switch purpose
        case 'collocation'
            if ~isdae
                options.type='gauss';
            else
                options.type='cheb';
            end                
        otherwise
            if ~isdae && degree<=options.maxdegree_linear && ~options.force_smooth
                options.type='equidistant';
            else
                options.type='cheb';
            end
    end
end
switch options.type
    case {'gauss','legendre'}
        grid=legpts(degree,[0,1]);
        grid=grid(:).';
    case {'cheb','chebyshev'}
        grid=dde_coll_chebxwt(degree,1:degree+1,[0,1]);
    case {'force_smooth','chebm1','chebyshev_smooth'}
        switch purpose
            case 'storage'
                grid=dde_coll_chebxwt(degree,1:degree+1,[0,1]);
            case 'collocation'        
                degree=degree-1;
                if isdae
                    grid=dde_coll_chebxwt(degree+1,2:degree+2,[0,1]);
                else
                    grid=dde_coll_chebxwt(degree,1:degree+1,[0,1]);
                end
        end
    case {'linear','equidistant'}
        grid=linspace(0,1,degree+1);
end
switch purpose
    case 'collocation'
        if length(grid)>degree
            grid=options.subset(grid);
        end
end
end
