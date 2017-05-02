function G=wf_enh(X,N,delta)

%WF_ENH: Wiener filter enhancement gain estimation
%    usage: G=WF_ENH(X,N,delta);

G = max((X-N)./X,delta);
end
