function supertitle(mytitle,varargin)
% When using subtitle('MY TITLE','PorpertyName','PropertyValue'...), or
% subtitle('MY TITLE') after a group of subplots, then it provides a title
% MY TITLE with any property used that is defined in the original title
% function in Matlab, but without affecting the titles, xlables and ylabels
% of any of the subplots. 
%
% Make sure use the function after the group of subplots.
%
% Example
% x = 0:0.01:6;
% subplot(221), plot(x,sin(x)), xlabel('x'), ylabel('sin(x)'), title('sin(x)')
% subplot(222), plot(x,cos(x)), xlabel('x'), ylabel('cos(x)'), title('cos(x)')
% subplot(223), plot(x,sin(2*x)), xlabel('x'), ylabel('sin(2x)'), title('sin(2x)')
% subplot(224), plot(x,cos(2*x)), xlabel('x'), ylabel('cos(2x)'), title('cos(2x)')
% subtitle('Single title on top','FontSize',12,'Color','r')
%
% Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)

axes('Units','Normal');
h = title(mytitle,varargin{1:length(varargin)});
set(gca,'visible','off')
set(h,'visible','on')
end

% Copyright (c) 2016, Shoaibur
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.